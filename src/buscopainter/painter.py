#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BuscoPainter T7.1-fixed — paint BUSCO hits onto an assembly and write GenBank.

Derived from buscoPainter.T7.1.py. Every change is listed in CHANGELOG_T7_FIXES.md.

What it does, per sample:
- Finds the BUSCO run_* directory and reads full_table.tsv (BUSCO v5/v6 layout).
- Reads the per-BUSCO GFFs in run_*/busco_sequences/{single_copy,multi_copy,fragmented}_busco_sequences/.
- Matches each non-Missing full_table row to exactly ONE gene model (GFF mRNA + its CDS rows)
  by BUSCO id, contig and coordinate overlap. Duplicated copies stay separate features.
- Writes one `gene` (mRNA span) and one `CDS` (joined CDS rows, original coordinates) per hit.
  When BUSCO dropped terminal miniprot segments (its full_table span is a run of whole CDS rows
  and the .faa protein covers only those), only BUSCO's segments are written (T7.1-fixed.3).
  GFF phase is NOT used to move boundaries; the 5' segment's phase becomes /codon_start.
- Attaches /translation only when the CDS is frame-consistent, a whole number of codons and has
  no internal stop codon. Otherwise the reason is written to a /note and the issues table.
- Nothing is clamped or silently dropped: every rejected or flagged hit is written to
  <prefix>_<sample>.annotation_issues.tsv (and rejected ones also to .skipped.tsv).
- --dry-run performs all parsing, matching and checks and writes the TSV reports, but no GenBank.
- Only BUSCO marker genes are painted; every record's source feature says so (T7.1-fixed.2).
- --per-contig also writes one single-record GenBank per contig: Artemis opens only the first
  record of a multi-record GenBank file (T7.1-fixed.2).
- A failed count reconciliation (exported + rejected != hit rows) gives exit code 3.

Counting conventions (written to the summary TSV with explicit column names):
- *_orthologs : number of distinct BUSCO ids per status (the BUSCO short_summary convention;
                C = S + D).
- *_hit_rows  : number of full_table rows per status (a Duplicated BUSCO has >= 2 rows).
- exported_genes / exported_CDS : one per hit row successfully placed on the assembly.
- cds_segments : number of GFF CDS rows joined into those CDS features.
"""

import argparse, sys, json, gzip, logging, re, time, csv
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
import pandas as pd

# optional deps
try:
    import yaml  # type: ignore
except Exception:
    yaml = None
try:
    import requests  # type: ignore
except Exception:
    requests = None

from urllib.parse import quote
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation

__version__ = "T7.1-fixed.3"

# --------------------------- Defaults ----------------------------------------
STATUS_DIRS = {
    "Complete": "single_copy_busco_sequences",
    "Duplicated": "multi_copy_busco_sequences",
    "Fragmented": "fragmented_busco_sequences",
}
KNOWN_STATUSES = ("Complete", "Duplicated", "Fragmented", "Missing")
FULL_TABLE_COLUMNS = ["Busco id", "Status", "Sequence", "Gene Start", "Gene End", "Strand",
                      "Score", "Length", "OrthoDB url", "Description"]
# full_table spans differ from the GFF mRNA span by up to 3 bp (BUSCO reports a 0-based start
# and excludes the stop codon). Larger differences are reported as a warning.
SPAN_TOLERANCE_BP = 3
TRANSL_TABLE = 1  # standard code; BUSCO genes are nuclear

ORTHODB_API_BASE = "https://data.orthodb.org/{version}/orthologs?id={og}"
ORTHODB_API_VERSION_DEFAULT = "v12"
ORTHODB_PAUSE_SECONDS_DEFAULT = 1.0

LOG_FORMAT_CONSOLE = "[%(levelname)s] %(message)s"
LOG_FORMAT_FILE = "%(asctime)s [%(levelname)s] %(message)s"

ISSUE_COLUMNS = ["sample", "busco_id", "status", "sequence", "ft_start", "ft_end", "ft_strand",
                 "mrna_id", "severity", "code", "detail"]

# --------------------------- Logging -----------------------------------------
logger = logging.getLogger("buscoPainter")
logger.setLevel(logging.INFO)
logger.handlers.clear()
_console = logging.StreamHandler(sys.stdout)
_console.setLevel(logging.INFO)
_console.setFormatter(logging.Formatter(LOG_FORMAT_CONSOLE))
logger.addHandler(_console)


def setup_file_logger(outdir: Path, prefix: str, log_level: str = "INFO"):
    outdir.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(outdir / f"{prefix}.busco_painter.log", mode="w", encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter(LOG_FORMAT_FILE))
    logger.addHandler(fh)
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))


class PainterError(Exception):
    """Fatal input problem for a sample."""


# --------------------------- Path checks -------------------------------------
def _ensure_no_spaces(path_str: Optional[str], label: str):
    if path_str and re.search(r"\s", path_str):
        raise PainterError(f"{label} contains whitespace: {path_str}")


# --------------------------- Config loading ----------------------------------
def load_config(path: Optional[str]) -> Dict:
    if not path:
        return {}
    _ensure_no_spaces(path, "--config")
    p = Path(path)
    if not p.exists():
        raise PainterError(f"--config file not found: {path}")
    text = p.read_text(encoding="utf-8")
    if p.suffix.lower() in (".yaml", ".yml"):
        if yaml is None:
            raise PainterError("YAML config provided but PyYAML not installed.")
        return yaml.safe_load(text)
    return json.loads(text)


# --------------------------- BUSCO helpers -----------------------------------
def _full_table_in(dirpath: Path) -> Optional[Path]:
    for name in ("full_table.tsv", "full_table.tsv.gz"):
        p = dirpath / name
        if p.exists():
            return p
    return None


def locate_busco_run_dir(base_dir: Path) -> Tuple[Path, Path]:
    ft_here = _full_table_in(base_dir)
    if ft_here:
        return base_dir, ft_here
    candidates = []
    if base_dir.exists():
        for child in sorted(base_dir.iterdir()):
            if child.is_dir() and child.name.startswith("run_"):
                ft = _full_table_in(child)
                if ft:
                    candidates.append((child, ft))
    if len(candidates) == 1:
        run_dir, ft = candidates[0]
        logger.info(f"Detected BUSCO run directory: {run_dir}")
        return run_dir, ft
    raise PainterError(f"Could not uniquely identify a BUSCO run_* directory under {base_dir} "
                       f"(found {len(candidates)})")


def get_status_gff_dirs(run_dir: Path) -> Dict[str, Path]:
    """Map BUSCO status -> folder holding that status' GFFs (BUSCO v5/v6 and older layouts)."""
    out = {}
    for status, sub in STATUS_DIRS.items():
        for cand in (run_dir / "busco_sequences" / sub, run_dir / sub,
                     run_dir / "busco_sequences" / sub.replace("_copy", "")):
            if cand.is_dir():
                out[status] = cand
                break
    return out


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def read_full_table(full_table_path: Path) -> pd.DataFrame:
    """Read BUSCO full_table.tsv. The header is a '# Busco id ...' comment line in v5/v6.
    Missing rows carry only id+status. Returns string columns plus integer Gene Start/End."""
    rows = []
    with _open_text(full_table_path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) > len(FULL_TABLE_COLUMNS):
                raise PainterError(f"{full_table_path}:{lineno}: {len(parts)} columns, "
                                   f"expected <= {len(FULL_TABLE_COLUMNS)}")
            parts += [""] * (len(FULL_TABLE_COLUMNS) - len(parts))
            if parts[1] not in KNOWN_STATUSES:
                raise PainterError(f"{full_table_path}:{lineno}: unknown status {parts[1]!r}")
            rows.append(parts + [lineno])
    df = pd.DataFrame(rows, columns=FULL_TABLE_COLUMNS + ["ft_line"], dtype=object)
    for c in ("Gene Start", "Gene End"):
        df[c] = pd.to_numeric(df[c].replace("", None), errors="coerce").astype("Int64")
    bad = df[(df["Status"] != "Missing") &
             (df["Gene Start"].isna() | df["Gene End"].isna() | (df["Sequence"] == "") |
              ~df["Strand"].isin(["+", "-"]))]
    if len(bad):
        raise PainterError(f"{full_table_path}: {len(bad)} non-Missing rows lack sequence/"
                           f"coordinates/strand (first at line {bad.iloc[0]['ft_line']})")
    return df


def busco_counts(df: pd.DataFrame) -> Dict[str, int]:
    """Ortholog-level counts (distinct BUSCO ids, BUSCO convention) and hit-row counts."""
    out = {}
    for st, key in (("Complete", "S"), ("Duplicated", "D"), ("Fragmented", "F"), ("Missing", "M")):
        sub = df[df["Status"] == st]
        out[f"busco_{key}_orthologs"] = int(sub["Busco id"].nunique())
        if st != "Missing":
            out[f"busco_{key}_hit_rows"] = int(len(sub))
    out["busco_C_orthologs"] = out["busco_S_orthologs"] + out["busco_D_orthologs"]
    out["busco_n_orthologs"] = int(df["Busco id"].nunique())
    return out


def read_short_summary(run_dir: Path) -> Optional[Dict[str, int]]:
    p = run_dir / "short_summary.json"
    if not p.exists():
        return None
    try:
        res = json.loads(p.read_text(encoding="utf-8"))["results"]
        return {"S": int(res["Single copy BUSCOs"]), "D": int(res["Multi copy BUSCOs"]),
                "F": int(res["Fragmented BUSCOs"]), "M": int(res["Missing BUSCOs"]),
                "n": int(res["n_markers"])}
    except Exception as e:  # noqa: BLE001
        logger.warning(f"Could not read {p}: {e}")
        return None


# --------------------------- OrthoDB optional -------------------------------
def fetch_orthodb_products_via_api(og_ids: Set[str], outdir: Path,
                                   api_version: str = ORTHODB_API_VERSION_DEFAULT,
                                   pause_seconds: float = ORTHODB_PAUSE_SECONDS_DEFAULT) -> Dict[str, str]:
    if requests is None:
        logger.warning("requests not installed; cannot fetch OrthoDB.")
        return {}
    cache_path = outdir / "orthodb_cache.json"
    cache = {}
    if cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text(encoding="utf-8"))
        except Exception:
            cache = {}
    session = requests.Session()
    product_map = {}
    for og in sorted(og_ids):
        entry = cache.get(og)
        if not entry:
            url = ORTHODB_API_BASE.format(version=api_version, og=quote(og))
            try:
                resp = session.get(url, timeout=15)
            except Exception as e:  # noqa: BLE001
                logger.warning(f"OrthoDB request failed for {og}: {e}")
                continue
            if resp.status_code != 200:
                logger.warning(f"OrthoDB HTTP {resp.status_code} for {og}")
                continue
            entry = resp.json()
            cache[og] = entry
            time.sleep(pause_seconds)
        descs = [(g.get("description") or "").strip() for g in entry.get("data", {}).get("genes", [])]
        descs = [d for d in descs if d]
        if descs:
            product_map[og] = descs[0]
    try:
        cache_path.write_text(json.dumps(cache, indent=2), encoding="utf-8")
    except Exception:
        pass
    return product_map


# --------------------------- GFF parsing -------------------------------------
def _parse_attrs(col9: str) -> Dict[str, str]:
    out = {}
    for kv in col9.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_busco_gff(gff_path: Path) -> Tuple[Dict[str, Dict], List[str]]:
    """Parse one BUSCO/miniprot GFF into gene models keyed by mRNA ID.
    Coordinates are kept as in the file (1-based, inclusive). Returns (models, problems)."""
    models: Dict[str, Dict] = {}
    children = defaultdict(list)
    bad_children = defaultdict(list)  # Parent -> unparseable CDS/stop_codon rows
    problems: List[str] = []
    with open(gff_path, "r", encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9:
                problems.append(f"line {lineno}: {len(p)} columns")
                continue
            seqid, _src, ftype, start, end, _score, strand, phase, attrs = p
            a = _parse_attrs(attrs)
            try:
                start, end = int(start), int(end)
            except ValueError:
                problems.append(f"line {lineno}: non-integer coordinates")
                if ftype in ("CDS", "stop_codon") and a.get("Parent"):
                    bad_children[a["Parent"]].append(f"{ftype} line {lineno} unparseable")
                continue
            rec = {"seqid": seqid, "type": ftype, "start": start, "end": end, "strand": strand,
                   "phase": phase, "attrs": a, "line": lineno}
            if start < 1 or end < start or strand not in ("+", "-"):
                why = (f"invalid interval {start}-{end}" if start < 1 or end < start
                       else f"invalid strand {strand!r}")
                problems.append(f"line {lineno}: {why}")
                if ftype in ("CDS", "stop_codon") and a.get("Parent"):
                    bad_children[a["Parent"]].append(f"{ftype} line {lineno} {why}")
                continue
            if ftype == "mRNA":
                mid = a.get("ID")
                if not mid or mid in models:
                    problems.append(f"line {lineno}: mRNA without unique ID")
                    continue
                models[mid] = {**rec, "id": mid, "cds": [], "stop_codon": []}
            elif ftype in ("CDS", "stop_codon"):
                par = a.get("Parent")
                if not par:
                    problems.append(f"line {lineno}: {ftype} without Parent")
                    continue
                if ftype == "CDS" and phase not in ("0", "1", "2"):
                    problems.append(f"line {lineno}: CDS phase {phase!r} not in 0/1/2")
                    bad_children[par].append(f"CDS line {lineno} phase {phase!r} not in 0/1/2")
                    continue
                children[par].append(rec)
    # A model that lost a CDS row at parse time would be exported with a segment missing:
    # mark it invalid so it is rejected and reported, never written as a partial gene.
    for par, why in bad_children.items():
        if par in models:
            models[par].setdefault("errors", []).extend(why)
    for par, recs in children.items():
        if par not in models:
            problems.append(f"{len(recs)} child rows reference unknown mRNA {par}")
            continue
        m = models[par]
        for r in recs:
            if r["seqid"] != m["seqid"] or r["strand"] != m["strand"]:
                m.setdefault("errors", []).append(f"{r['type']} line {r['line']} seqid/strand differs from mRNA")
            elif r["start"] < m["start"] or r["end"] > m["end"]:
                m.setdefault("errors", []).append(f"{r['type']} line {r['line']} outside mRNA span")
            if r["type"] == "CDS":
                m["cds"].append(r)
            else:
                m["stop_codon"].append(r)
    for m in models.values():
        # transcription order: ascending on +, descending on -
        m["cds"].sort(key=lambda r: r["start"], reverse=(m["strand"] == "-"))
        asc = sorted(m["cds"], key=lambda r: r["start"])
        for a_, b_ in zip(asc, asc[1:]):
            if b_["start"] <= a_["end"]:
                m.setdefault("errors", []).append(f"overlapping CDS rows at lines {a_['line']},{b_['line']}")
        if not m["cds"]:
            m.setdefault("errors", []).append("mRNA has no CDS rows")
    return models, problems


def frame_consistency(cds_rows: List[Dict]) -> bool:
    """True when each segment's GFF phase equals the phase implied by the preceding segments
    (i.e. the joined CDS is one continuous reading frame). cds_rows in transcription order."""
    if not cds_rows:
        return False
    first_phase = int(cds_rows[0]["phase"])
    consumed = 0
    for i, r in enumerate(cds_rows):
        if i:
            expected = (3 - ((consumed - first_phase) % 3)) % 3
            if int(r["phase"]) != expected:
                return False
        consumed += r["end"] - r["start"] + 1
    return True


def build_cds_location(cds_rows: List[Dict]):
    """GFF 1-based inclusive -> Biopython 0-based half-open. No boundary changes."""
    strand = 1 if cds_rows[0]["strand"] == "+" else -1
    parts = [SimpleLocation(r["start"] - 1, r["end"], strand=strand) for r in cds_rows]
    return parts[0] if len(parts) == 1 else CompoundLocation(parts)


def evaluate_translation(contig_seq: Seq, cds_loc, codon_start: int, frame_ok: bool,
                         blocker: str = "") -> Tuple[Optional[str], str, Dict]:
    """Return (protein or None, reason, facts). Only a clean, single-frame translation is returned.
    `blocker` is a model-level reason (e.g. miniprot frameshift) that forbids a translation."""
    nt = cds_loc.extract(contig_seq)
    coding = nt[codon_start - 1:]
    facts = {"starts_with_ATG": str(nt[:3]).upper() == "ATG" and codon_start == 1,
             "ends_with_stop": False, "internal_stops": None}
    if blocker:
        return None, blocker, facts
    if not frame_ok:
        return None, "CDS segments are not in one continuous reading frame (GFF phases)", facts
    if len(coding) == 0:
        return None, "empty CDS", facts
    if len(coding) % 3:
        return None, f"coding length after codon_start ({len(coding)} nt) is not a multiple of 3", facts
    prot = str(coding.translate(table=TRANSL_TABLE))
    if prot.endswith("*"):
        facts["ends_with_stop"] = True
        prot = prot[:-1]
    n_int = prot.count("*")
    facts["internal_stops"] = n_int
    if n_int:
        return None, f"{n_int} internal stop codon(s) in the GFF-defined frame", facts
    if not prot:
        return None, "empty translation", facts
    if set(prot) == {"X"}:
        return None, "translation is entirely ambiguous (X)", facts
    return prot, "", facts


# --------------------------- Feature writing ---------------------------------
def scope_note(lineage: str) -> str:
    """Stated on every record so a viewer does not read unannotated ORFs as missed BUSCOs."""
    return (f"BuscoPainter {__version__}: only BUSCO {lineage} marker genes are annotated. "
            "This is not a whole-genome gene annotation; unannotated ORFs are expected.")


def add_source_feature(rec, organism: str, isolate: Optional[str], topology: str,
                       note: Optional[str] = None):
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = topology
    if organism:
        rec.annotations["organism"] = organism
    if not any(f.type == "source" for f in rec.features):
        qualifiers = {"organism": [organism] if organism else [], "mol_type": ["genomic DNA"]}
        if isolate:
            qualifiers["isolate"] = [isolate]
        if note:
            qualifiers["note"] = [note]
        src = SeqFeature(SimpleLocation(0, len(rec), strand=1), type="source", qualifiers=qualifiers)
        rec.features.insert(0, src)


def split_cds_by_full_table_span(cds_rows: List[Dict], lo: int, hi: int):
    """BUSCO (miniprot mode) drops terminal miniprot segments it does not accept: full_table and
    the .faa protein then cover only the remaining segments, while the GFF keeps the whole model.
    Returns (kept, dropped) when the full_table span lo..hi is exactly a run of whole CDS rows
    (within SPAN_TOLERANCE_BP) and every other row lies wholly outside it; otherwise None."""
    tol = SPAN_TOLERANCE_BP
    kept = [r for r in cds_rows if r["start"] >= lo - tol and r["end"] <= hi + tol]
    dropped = [r for r in cds_rows if r["end"] < lo or r["start"] > hi]
    if not kept or len(kept) + len(dropped) != len(cds_rows):
        return None
    if abs(min(r["start"] for r in kept) - lo) > tol or abs(max(r["end"] for r in kept) - hi) > tol:
        return None
    return kept, dropped


def _loc_string(rows: List[Dict], strand: str) -> str:
    asc = sorted(rows, key=lambda r: r["start"])
    s = ",".join(f"{r['start']}..{r['end']}" for r in asc)
    s = f"join({s})" if len(asc) > 1 else s
    return f"complement({s})" if strand == "-" else s


def match_rows_to_models(rows: List[Dict], models: Dict[str, Dict]) -> Tuple[Dict[int, str], Dict[int, str]]:
    """One-to-one assignment of full_table rows to GFF mRNAs of the same BUSCO id.
    Candidates: same contig and overlapping span. Best = smallest |start diff|+|end diff|.
    Returns ({row_index: mrna_id}, {row_index: failure reason})."""
    pairs = []
    for i, r in enumerate(rows):
        lo, hi = sorted((int(r["Gene Start"]), int(r["Gene End"])))
        for mid, m in models.items():
            if m["seqid"] == r["Sequence"] and m["start"] <= hi and lo <= m["end"]:
                pairs.append((abs(m["start"] - lo) + abs(m["end"] - hi), i, mid))
    pairs.sort()
    assigned, used, failed = {}, set(), {}
    cand_by_row = defaultdict(list)
    for d, i, mid in pairs:
        cand_by_row[i].append((d, mid))
    for d, i, mid in pairs:
        if i in assigned or mid in used:
            continue
        free_equal = [m2 for d2, m2 in cand_by_row[i] if d2 == d and m2 not in used]
        if len(free_equal) > 1:
            failed[i] = f"AMBIGUOUS_MATCH:{','.join(sorted(free_equal))}"
            continue
        assigned[i] = mid
        used.add(mid)
    for i in range(len(rows)):
        if i not in assigned and i not in failed:
            failed[i] = ("MODEL_ALREADY_ASSIGNED" if cand_by_row[i] else "NO_OVERLAPPING_GFF_MODEL")
    return assigned, failed


def add_busco_features(records: Dict[str, "SeqRecord"],
                       busco_rows: pd.DataFrame,
                       status_dirs: Dict[str, Path],
                       run_dir: Path,
                       orthodb_map: Dict[str, str],
                       sample_name: str,
                       lineage: str,
                       write_features: bool = True) -> Tuple[List[Dict], List[Dict]]:
    """Attach gene+CDS per BUSCO hit. Returns (feature_table_rows, issues)."""
    issues: List[Dict] = []
    feats_out: List[Dict] = []
    locus_tags: Set[str] = set()

    def issue(r, severity, code, detail="", mrna_id=""):
        issues.append({"sample": sample_name, "busco_id": r["Busco id"], "status": r["Status"],
                       "sequence": r.get("Sequence", ""), "ft_start": r.get("Gene Start", ""),
                       "ft_end": r.get("Gene End", ""), "ft_strand": r.get("Strand", ""),
                       "mrna_id": mrna_id, "severity": severity, "code": code, "detail": detail})

    present = busco_rows[busco_rows["Status"] != "Missing"]
    for busco_id, grp in present.groupby("Busco id", sort=False):
        rows = grp.to_dict("records")
        statuses = {r["Status"] for r in rows}
        if len(statuses) > 1:
            for r in rows:
                issue(r, "REJECTED", "CONFLICTING_STATUS", f"statuses {sorted(statuses)} for one BUSCO id")
            continue
        status = rows[0]["Status"]
        if status == "Complete" and len(rows) > 1:
            for r in rows:
                issue(r, "WARNING", "MULTIPLE_ROWS_FOR_SINGLE_COPY", f"{len(rows)} rows")

        # Locate GFF: expected folder for the status first, then any other folder (reported).
        gff = None
        exp_dir = status_dirs.get(status)
        if exp_dir and (exp_dir / f"{busco_id}.gff").exists():
            gff = exp_dir / f"{busco_id}.gff"
        else:
            for st, d in status_dirs.items():
                if (d / f"{busco_id}.gff").exists():
                    gff = d / f"{busco_id}.gff"
                    for r in rows:
                        issue(r, "WARNING", "GFF_IN_UNEXPECTED_FOLDER", str(gff.relative_to(run_dir)))
                    break
        if gff is None:
            for r in rows:
                issue(r, "REJECTED", "MISSING_GFF", f"{busco_id}.gff not found")
            continue

        models, problems = parse_busco_gff(gff)
        for pr in problems:
            issue(rows[0], "WARNING", "GFF_PARSE_PROBLEM", f"{gff.name}: {pr}")

        assigned, failed = match_rows_to_models(rows, models)
        for i, reason in failed.items():
            issue(rows[i], "REJECTED", reason.split(":")[0], f"{gff.name} {reason}")
        unused = set(models) - set(assigned.values())
        for mid in sorted(unused):
            m = models[mid]
            issue(rows[0], "INFO", "GFF_MODEL_NOT_IN_FULL_TABLE",
                  f"{gff.name} {mid} {m['seqid']}:{m['start']}-{m['end']}({m['strand']}) not exported", mid)

        gff_rel = str(gff.relative_to(run_dir))
        for i, mid in sorted(assigned.items()):
            r, m = rows[i], models[mid]
            contig = r["Sequence"]
            if m.get("errors"):
                issue(r, "REJECTED", "INVALID_GFF_MODEL", "; ".join(m["errors"]), mid)
                continue
            if m["strand"] != r["Strand"]:
                issue(r, "REJECTED", "STRAND_MISMATCH",
                      f"full_table {r['Strand']} vs GFF {m['strand']}", mid)
                continue
            if contig not in records:
                issue(r, "REJECTED", "MISSING_CONTIG_IN_FASTA", contig, mid)
                continue
            rec = records[contig]
            seq_len = len(rec)
            if m["end"] > seq_len or any(c["end"] > seq_len for c in m["cds"]):
                issue(r, "REJECTED", "OUT_OF_RANGE", f"{m['start']}-{m['end']} > contig length {seq_len}", mid)
                continue

            cds_rows = m["cds"]
            gene_start, gene_end = m["start"], m["end"]
            flags = []
            blockers = []
            dropped = []
            lo, hi = sorted((int(r["Gene Start"]), int(r["Gene End"])))
            if abs(m["start"] - lo) > SPAN_TOLERANCE_BP or abs(m["end"] - hi) > SPAN_TOLERANCE_BP:
                split = split_cds_by_full_table_span(cds_rows, lo, hi)
                if split:
                    # T7.1-fixed.3: export BUSCO's prediction (full_table + .faa), not the
                    # untrimmed miniprot model.
                    cds_rows, dropped = split
                    gene_start = min(c["start"] for c in cds_rows)
                    gene_end = max(c["end"] for c in cds_rows)
                    flags.append("FULLTABLE_SPAN_TRIMMED")
                    issue(r, "WARNING", "FULLTABLE_SPAN_TRIMMED",
                          f"full_table {lo}-{hi} vs GFF mRNA {m['start']}-{m['end']}; exported the "
                          f"{len(cds_rows)} CDS segment(s) inside the BUSCO span; not exported: "
                          + ",".join(f"{c['start']}-{c['end']}" for c in dropped), mid)
                else:
                    flags.append("FULLTABLE_SPAN_DIFFERS_FROM_GFF")
                    issue(r, "WARNING", "FULLTABLE_SPAN_DIFFERS_FROM_GFF",
                          f"full_table {lo}-{hi} vs GFF mRNA {m['start']}-{m['end']}; GFF span exported", mid)
                    blockers.append("BUSCO full_table span covers only part of the GFF model")
            codon_start = int(cds_rows[0]["phase"]) + 1
            frame_ok = frame_consistency(cds_rows)
            cds_loc = build_cds_location(cds_rows)
            strand_num = 1 if m["strand"] == "+" else -1
            gene_loc = SimpleLocation(gene_start - 1, gene_end, strand=strand_num)
            # miniprot tags Frameshift/StopCodon on the affected CDS rows (and on the mRNA); for a
            # trimmed model only the exported rows count.
            model_attrs = [] if dropped else [m["attrs"]]
            if any("Frameshift" in a for a in model_attrs + [c["attrs"] for c in cds_rows]):
                blockers.append("miniprot reports frameshift(s) inside the alignment")
            if any("StopCodon" in a for a in model_attrs + [c["attrs"] for c in cds_rows]):
                blockers.append("miniprot reports in-frame stop codon(s)")
            prot, why_not, facts = evaluate_translation(rec.seq, cds_loc, codon_start, frame_ok,
                                                        "; ".join(blockers))
            gaps = []
            asc = sorted(cds_rows, key=lambda c: c["start"])
            for a_, b_ in zip(asc, asc[1:]):
                gaps.append(b_["start"] - a_["end"] - 1)
            if len(cds_rows) > 1:
                flags.append("MULTI_SEGMENT_MODEL")
                issue(r, "WARNING", "MULTI_SEGMENT_MODEL",
                      f"{len(cds_rows)} CDS segments, gaps {gaps} bp (unverified spliced alignment)", mid)
            if not frame_ok:
                flags.append("FRAME_INCONSISTENT")
                issue(r, "WARNING", "FRAME_INCONSISTENT",
                      "phases " + ",".join(c["phase"] for c in cds_rows), mid)
            if "Frameshift" in m["attrs"]:
                flags.append("MINIPROT_FRAMESHIFT")
            if "StopCodon" in m["attrs"]:
                flags.append("MINIPROT_STOPCODON")
            if prot is None:
                flags.append("NO_TRANSLATION")
                issue(r, "WARNING", "TRANSLATION_OMITTED", why_not, mid)

            locus_tag = f"{sample_name}_{busco_id}_{mid}"
            if locus_tag in locus_tags:
                issue(r, "REJECTED", "DUPLICATE_LOCUS_TAG", locus_tag, mid)
                continue
            locus_tags.add(locus_tag)

            desc = str(r.get("Description") or "").strip()
            product = desc or orthodb_map.get(busco_id) or "hypothetical protein"
            ev = (f"BUSCO {busco_id} ({status}), lineage {lineage}; full_table line {r['ft_line']}: "
                  f"{contig}:{r['Gene Start']}-{r['Gene End']}({r['Strand']}) score={r['Score']} "
                  f"length={r['Length']}")
            gffev = (f"GFF evidence: {gff_rel} mRNA {mid} {m['seqid']}:{m['start']}-{m['end']}({m['strand']}); "
                     f"Identity={m['attrs'].get('Identity', 'NA')}; Target={m['attrs'].get('Target', 'NA')}")
            notes = [ev, gffev]
            if dropped:
                notes.append(f"BUSCO's prediction (full_table span and protein) covers {len(cds_rows)} of the "
                             f"{len(m['cds'])} miniprot CDS segments; segments outside it were not exported: "
                             + ",".join(f"{c['start']}..{c['end']}" for c in dropped))
            if len(cds_rows) > 1:
                notes.append(f"CDS joins {len(cds_rows)} miniprot segments separated by {gaps} bp; "
                             "cis-introns are rare in trypanosomatids, so this structure is unverified")
            if not frame_ok:
                notes.append("GFF phases are not consistent with one reading frame across segments "
                             "(possible frameshift); exported as given")
            if prot is None:
                notes.append(f"translation omitted: {why_not}")
            gene_q = {"locus_tag": [locus_tag], "note": [f"BUSCO {busco_id} {status}; {gffev}"]}
            cds_q = {"locus_tag": [locus_tag], "product": [product], "codon_start": [str(codon_start)],
                     "transl_table": [str(TRANSL_TABLE)], "note": notes}
            if prot is not None:
                cds_q["translation"] = [prot]
            if write_features:
                rec.features.append(SeqFeature(gene_loc, type="gene", qualifiers=gene_q))
                rec.features.append(SeqFeature(cds_loc, type="CDS", qualifiers=cds_q))
            feats_out.append({
                "sample": sample_name, "locus_tag": locus_tag, "busco_id": busco_id, "status": status,
                "contig": contig, "contig_length": seq_len, "strand": m["strand"], "mrna_id": mid,
                "gene_start": gene_start, "gene_end": gene_end,
                "cds_location": _loc_string(cds_rows, m["strand"]), "cds_segments": len(cds_rows),
                "codon_start": codon_start, "frame_consistent": frame_ok,
                "translation_attached": prot is not None, "translation_length": len(prot) if prot else 0,
                "translation_issue": why_not, "starts_with_ATG": facts["starts_with_ATG"],
                "ends_with_stop": facts["ends_with_stop"], "flags": ";".join(flags),
                "ft_line": r["ft_line"], "ft_start": r["Gene Start"], "ft_end": r["Gene End"],
                "ft_score": r["Score"], "ft_length": r["Length"], "gff_file": gff_rel,
                "product": product,
            })
    return feats_out, issues


def _sort_features(rec):
    src = [f for f in rec.features if f.type == "source"]
    rest = [f for f in rec.features if f.type != "source"]
    order = {"gene": 0, "CDS": 1}
    rest.sort(key=lambda f: (int(f.location.start), -int(f.location.end), order.get(f.type, 9),
                             f.qualifiers.get("locus_tag", [""])[0]))
    rec.features = src + rest


def write_genbank(out_path: Path, records: List["SeqRecord"]):
    with open(out_path, "w", encoding="utf-8") as out:
        SeqIO.write(records, out, "genbank")
    logger.info(f"Wrote GenBank: {out_path}")


def safe_file_part(name: str) -> str:
    """Contig id -> file-name-safe string (path separators and other odd characters become '_')."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", name).lstrip(".") or "_"


def write_per_contig_genbanks(outdir: Path, prefix: str, sample: str, records: Dict[str, "SeqRecord"]):
    """One single-record GenBank per contig. Artemis opens only the first record of a
    multi-record GenBank file, so these are the files to open there."""
    d = outdir / f"{prefix}_{sample}_per_contig"
    d.mkdir(parents=True, exist_ok=True)
    seen: Dict[str, str] = {}
    for cid, rec in records.items():
        fn = f"{prefix}_{sample}_contig-{safe_file_part(cid)}.gbk"
        if fn in seen:
            raise PainterError(f"Contig ids {seen[fn]!r} and {cid!r} map to the same file name {fn}")
        seen[fn] = cid
        with open(d / fn, "w", encoding="utf-8") as fh:
            SeqIO.write([rec], fh, "genbank")
    logger.info(f"Wrote {len(records)} per-contig GenBank files to {d}")


# --------------------------- Pre-scan + ID gather ----------------------------
def prescan_samples_and_collect_ids(samples: List[Dict], global_topology: str,
                                    global_organism: Optional[str]) -> Tuple[List[Dict], Set[str]]:
    all_ids: Set[str] = set()
    updated: List[Dict] = []
    for s in samples:
        name = s.get("name") or Path(s["fasta"]).stem
        fasta_path = Path(s["fasta"])
        busco_base = Path(s["busco_dir"])
        _ensure_no_spaces(str(fasta_path), f"sample[{name}] fasta")
        _ensure_no_spaces(str(busco_base), f"sample[{name}] busco_dir")
        if not fasta_path.exists():
            raise PainterError(f"FASTA not found: {fasta_path}")
        if not busco_base.exists():
            raise PainterError(f"BUSCO directory does not exist: {busco_base}")
        run_dir, full_table = locate_busco_run_dir(busco_base)
        df = read_full_table(full_table)
        all_ids.update(df.loc[df["Status"] != "Missing", "Busco id"].astype(str))
        lineage = s.get("lineage") or (run_dir.name[4:] if run_dir.name.startswith("run_") else "NA")
        updated.append({**s, "name": name, "run_dir": run_dir, "full_table": full_table,
                        "full_table_df": df, "lineage": lineage,
                        "topology": s.get("topology", global_topology),
                        "organism": s.get("organism", global_organism)})
    return updated, all_ids


# --------------------------- Pipeline ----------------------------------------
def load_fasta(fasta: Path) -> Dict[str, "SeqRecord"]:
    records = {}
    for r in SeqIO.parse(str(fasta), "fasta"):
        if r.id in records:
            raise PainterError(f"Duplicate FASTA id {r.id} in {fasta}")
        records[r.id] = r
    if not records:
        raise PainterError(f"No sequences in {fasta}")
    return records


def process_sample(sample: Dict, outdir: Path, prefix: str, orthodb: Dict[str, str],
                   dry_run: bool = False) -> Dict:
    name, fasta, run_dir, df = sample["name"], Path(sample["fasta"]), sample["run_dir"], sample["full_table_df"]
    logger.info(f"Processing {name} ({'DRY RUN' if dry_run else 'write'}) run_dir={run_dir}")
    status_dirs = get_status_gff_dirs(run_dir)
    for st in ("Complete", "Duplicated", "Fragmented"):
        n = int((df["Status"] == st).sum())
        if n and st not in status_dirs:
            logger.warning(f"{name}: {n} {st} rows but no {STATUS_DIRS[st]} folder")
    logger.info(f"{name}: GFF folders: " + ", ".join(f"{k}={v}" for k, v in status_dirs.items()))

    records = load_fasta(fasta)
    for r in records.values():
        add_source_feature(r, sample.get("organism") or "Unknown", sample.get("isolate"),
                           sample.get("topology", "linear"), note=scope_note(sample["lineage"]))

    feats, issues = add_busco_features(records, df, status_dirs, run_dir, orthodb, name,
                                       sample["lineage"], write_features=not dry_run)
    gbk_prefix = sample.get("prefix") or prefix

    rejected = [i for i in issues if i["severity"] == "REJECTED"]
    for i in rejected:
        logger.warning(f"{name}: REJECTED {i['busco_id']} {i['status']} {i['sequence']}: {i['code']} {i['detail']}")
    pd.DataFrame(issues, columns=ISSUE_COLUMNS).to_csv(
        outdir / f"{gbk_prefix}_{name}.annotation_issues.tsv", sep="\t", index=False)
    pd.DataFrame(rejected, columns=ISSUE_COLUMNS).to_csv(
        outdir / f"{gbk_prefix}_{name}.skipped.tsv", sep="\t", index=False)
    feat_df = pd.DataFrame(feats)
    feat_df.to_csv(outdir / f"{gbk_prefix}_{name}.features.tsv", sep="\t", index=False)

    if not dry_run:
        for rec in records.values():
            _sort_features(rec)
        write_genbank(outdir / f"{gbk_prefix}_{name}.gbk", list(records.values()))
        if sample.get("export_per_contig", False):
            write_per_contig_genbanks(outdir, gbk_prefix, name, records)
        elif len(records) > 1:
            logger.info(f"{name}: the GenBank file holds {len(records)} records (one per contig). "
                        "Artemis shows only the first record; use --per-contig for Artemis.")
    else:
        logger.info(f"{name}: DRY RUN — no GenBank written")

    # Per-contig summary: hit rows per status located on each contig + exported features.
    present = df[df["Status"] != "Missing"]
    hit = present.groupby(["Sequence", "Status"]).size().unstack(fill_value=0)
    contig_stats = []
    for cid, rec in records.items():
        fc = feat_df[feat_df["contig"] == cid] if len(feat_df) else feat_df
        contig_stats.append({
            "contig": cid, "length": len(rec),
            **{f"{st}_hit_rows": int(hit.loc[cid].get(st, 0)) if cid in hit.index else 0
               for st in ("Complete", "Duplicated", "Fragmented")},
            "exported_genes": int(len(fc)), "exported_CDS": int(len(fc)),
            "cds_segments": int(fc["cds_segments"].sum()) if len(fc) else 0,
        })
    cs = pd.DataFrame(contig_stats)
    total = {"contig": "TOTAL", **{c: int(cs[c].sum()) for c in cs.columns if c != "contig"}}
    cs = pd.concat([cs, pd.DataFrame([total])], ignore_index=True)
    cs.to_csv(outdir / f"{gbk_prefix}_{name}.contig_summary.tsv", sep="\t", index=False)

    counts = busco_counts(df)
    ss = read_short_summary(run_dir)
    if ss:
        agree = (ss["S"], ss["D"], ss["F"], ss["M"], ss["n"]) == (
            counts["busco_S_orthologs"], counts["busco_D_orthologs"], counts["busco_F_orthologs"],
            counts["busco_M_orthologs"], counts["busco_n_orthologs"])
        if not agree:
            logger.warning(f"{name}: full_table ortholog counts disagree with short_summary.json {ss}")
    else:
        agree = None
    summary = {
        "sample": name, "fasta_contigs": len(records), **counts,
        "short_summary_agrees": agree,
        "exported_genes": int(len(feat_df)), "exported_CDS": int(len(feat_df)),
        "cds_segments": int(feat_df["cds_segments"].sum()) if len(feat_df) else 0,
        "multi_segment_models": int((feat_df["cds_segments"] > 1).sum()) if len(feat_df) else 0,
        "translations_attached": int(feat_df["translation_attached"].sum()) if len(feat_df) else 0,
        "rejected_hits": len(rejected),
        "warnings": sum(1 for i in issues if i["severity"] == "WARNING"),
        "dry_run": dry_run,
    }
    hit_rows = counts["busco_S_hit_rows"] + counts["busco_D_hit_rows"] + counts["busco_F_hit_rows"]
    summary["reconciled"] = summary["exported_genes"] + len(rejected) == hit_rows
    if not summary["reconciled"]:
        logger.error(f"{name}: exported ({summary['exported_genes']}) + rejected ({len(rejected)}) "
                     f"!= non-Missing hit rows ({hit_rows})")
    logger.info(f"{name}: hit rows={hit_rows} exported={summary['exported_genes']} "
                f"rejected={len(rejected)} warnings={summary['warnings']}")
    return summary


# --------------------------- Main --------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description=f"BuscoPainter {__version__}")
    p.add_argument("--config", required=True)
    p.add_argument("--outdir")
    p.add_argument("--prefix")
    p.add_argument("--samples", nargs="+", help="only process these sample names")
    p.add_argument("--dry-run", action="store_true",
                   help="parse, match and validate; write TSV reports but no GenBank")
    p.add_argument("--per-contig", action="store_true",
                   help="also write one GenBank file per contig (Artemis opens only the first "
                        "record of a multi-record file); same as export_per_contig in the config")
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args(argv)
    try:
        cfg = load_config(args.config)
        outdir = Path(args.outdir or cfg.get("outdir", "results"))
        outdir.mkdir(parents=True, exist_ok=True)
        prefix = args.prefix or cfg.get("prefix", "Project")
        setup_file_logger(outdir, prefix, args.log_level)
        logger.info(f"BuscoPainter {__version__}; config={args.config}; outdir={outdir}; dry_run={args.dry_run}")
        samples_cfg = cfg.get("samples") or []
        if args.samples:
            samples_cfg = [s for s in samples_cfg if s.get("name") in set(args.samples)]
        if not samples_cfg:
            raise PainterError("No samples to process.")
        samples = [{**s, "export_per_contig": args.per_contig or
                    s.get("export_per_contig", cfg.get("export_per_contig", False))}
                   for s in samples_cfg]
        prepared, ids = prescan_samples_and_collect_ids(samples, cfg.get("topology", "linear"), cfg.get("organism"))
        orthodb_map = fetch_orthodb_products_via_api(ids, outdir) if cfg.get("use_orthodb_api", False) else {}
        summary = [process_sample(s, outdir, prefix, orthodb_map, dry_run=args.dry_run) for s in prepared]
        pd.DataFrame(summary).to_csv(outdir / f"{prefix}.busco_summary.tsv", sep="\t", index=False)
        if not all(s["reconciled"] for s in summary):
            logger.error("Count reconciliation failed for at least one sample (see above); exit code 3.")
            return 3
        logger.info("Done.")
        return 0
    except PainterError as e:
        logger.error(str(e))
        return 2


if __name__ == "__main__":
    sys.exit(main())
