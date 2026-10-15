#!/usr/bin/env python3
"""
Whole-assembly Artemis export: indexed FASTA + indexed GFF3.

Artemis opens only the first record of a multi-record GenBank file. It can, however, open a
multi-FASTA that has a samtools ``.fai`` index together with a BGZF-compressed, tabix-indexed
GFF3. It then shows a contig drop-down in the entry toolbar, and the GFF3 features follow the
selected contig. This module writes that pair from a painter GenBank file:

  <outdir>/<prefix>_<sample>_artemis/
      <fasta file name>              byte-identical copy of the assembly FASTA
      <fasta file name>.fai          samtools-compatible FASTA index
      <prefix>_<sample>.busco.gff3.gz      BGZF-compressed GFF3 (all contigs)
      <prefix>_<sample>.busco.gff3.gz.tbi  tabix index (needs the `tabix` program)

The GFF3 is converted from the GenBank written by ``buscopainter``. It does not read BUSCO again
and does not change any gene model: every gene/CDS pair (same locus_tag) becomes one gene, one
mRNA and one CDS row per CDS segment, with the same coordinates, strand and joins. The sequences
are never concatenated or edited; the FASTA is copied byte for byte and checked against the
GenBank sequences.

Usage:
  buscopainter-artemis --config config.json --outdir results            # all samples
  buscopainter-artemis --config config.json --outdir results --export-root artemis_view
  buscopainter-artemis --gbk results/p_S1.gbk --fasta S1.fasta --outdir results/artemis
"""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import quote

from Bio import SeqIO, bgzf

from buscopainter import __version__

SOURCE = "BuscoPainter"
# GFF3 column 9: escape the reserved characters (and anything non-printable) in values.
_SAFE = " !\"#$'()*+-./0123456789:<>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
_BUSCO_NOTE = re.compile(r"^BUSCO (\S+) \((\w+)\), lineage (\S+);")


class ExportError(Exception):
    pass


def gff_escape(value: str) -> str:
    return quote(str(value), safe=_SAFE)


def _attrs(pairs: List[Tuple[str, object]]) -> str:
    out = []
    for k, v in pairs:
        vals = v if isinstance(v, list) else [v]
        vals = [gff_escape(x) for x in vals if x is not None and str(x) != ""]
        if vals:
            out.append(f"{k}={','.join(vals)}")
    return ";".join(out)


# ------------------------------------------------------------------ FASTA ---
def fasta_index(path: Path) -> List[Tuple[str, int, int, int, int]]:
    """samtools faidx rules: per record (name, length, offset, bases/line, bytes/line).
    Every line of a record except the last must have the same length."""
    entries = []
    name = None
    with open(path, "rb") as fh:
        offset = 0
        cur = None
        for raw in fh:
            line_bytes = len(raw)
            if raw.startswith(b">"):
                if cur:
                    entries.append(_close(cur, path))
                name = raw[1:].split()[0].decode() if raw[1:].split() else ""
                if not name:
                    raise ExportError(f"{path}: empty FASTA header at byte {offset}")
                cur = {"name": name, "len": 0, "off": offset + line_bytes, "lb": None, "lw": None,
                       "short_seen": False}
            elif cur is not None:
                seq = raw.rstrip(b"\r\n")
                if not seq:
                    cur["short_seen"] = True  # a blank line may only end the record
                    offset += line_bytes
                    continue
                if cur["short_seen"]:
                    raise ExportError(f"{path}: record {cur['name']} has lines of different "
                                      "lengths; samtools faidx cannot index it")
                if cur["lb"] is None:
                    cur["lb"], cur["lw"] = len(seq), line_bytes
                elif len(seq) != cur["lb"] or line_bytes != cur["lw"]:
                    if len(seq) > cur["lb"]:
                        raise ExportError(f"{path}: record {cur['name']} has lines of different "
                                          "lengths; samtools faidx cannot index it")
                    cur["short_seen"] = True
                cur["len"] += len(seq)
            offset += line_bytes
        if cur:
            entries.append(_close(cur, path))
    names = [e[0] for e in entries]
    if len(set(names)) != len(names):
        raise ExportError(f"{path}: duplicate FASTA ids")
    return entries


def _close(cur, path):
    if cur["lb"] is None:
        raise ExportError(f"{path}: record {cur['name']} has no sequence")
    return (cur["name"], cur["len"], cur["off"], cur["lb"], cur["lw"])


def write_fai(entries, out: Path):
    with open(out, "w") as fh:
        for e in entries:
            fh.write("\t".join(map(str, e)) + "\n")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------- GenBank ---
def _segments(loc) -> List[Tuple[int, int]]:
    """1-based inclusive (start, end) of each part, in transcription order."""
    return [(int(p.start) + 1, int(p.end)) for p in loc.parts]


def gff_rows_for_record(rec) -> Tuple[List[Tuple[int, int, str]], int]:
    """GFF3 lines for one GenBank record: [(start, type_rank, line)], number of gene models."""
    genes: Dict[str, object] = {}
    cdss: Dict[str, object] = {}
    for f in rec.features:
        if f.type not in ("gene", "CDS"):
            continue
        lt = f.qualifiers.get("locus_tag", [None])[0]
        if not lt:
            raise ExportError(f"{rec.id}: {f.type} at {f.location} has no locus_tag")
        target = genes if f.type == "gene" else cdss
        if lt in target:
            raise ExportError(f"{rec.id}: duplicate {f.type} locus_tag {lt}")
        target[lt] = f
    if set(genes) != set(cdss):
        raise ExportError(f"{rec.id}: gene/CDS locus_tags differ: "
                          f"{sorted(set(genes) ^ set(cdss))[:5]}")
    rows = []
    for lt, g in genes.items():
        c = cdss[lt]
        strand = {1: "+", -1: "-"}.get(g.location.strand)
        if strand is None or c.location.strand != g.location.strand:
            raise ExportError(f"{rec.id}: {lt} has no consistent strand")
        notes = c.qualifiers.get("note", [])
        m = _BUSCO_NOTE.match(notes[0]) if notes else None
        if not m:
            raise ExportError(f"{rec.id}: {lt} CDS note does not start with the BUSCO evidence line")
        busco_id, status, lineage = m.groups()
        gs, ge = int(g.location.start) + 1, int(g.location.end)
        segs = _segments(c.location)
        codon_start = int(c.qualifiers.get("codon_start", ["1"])[0])
        mrna = f"{lt}.mRNA"
        common = [("busco_id", busco_id), ("busco_status", status), ("lineage", lineage)]
        product = c.qualifiers.get("product", [""])[0]
        rows.append((gs, 0, "\t".join([rec.id, SOURCE, "gene", str(gs), str(ge), ".", strand, ".",
                     _attrs([("ID", lt), ("Name", busco_id), ("locus_tag", lt)] + common)])))
        rows.append((gs, 1, "\t".join([rec.id, SOURCE, "mRNA", str(gs), str(ge), ".", strand, ".",
                     _attrs([("ID", mrna), ("Parent", lt), ("Name", busco_id),
                             ("locus_tag", lt), ("product", product)])])))
        # GFF3 phase per segment, derived from /codon_start and the segment lengths in
        # transcription order (the reading frame the GenBank join defines).
        done = 0
        for i, (s, e) in enumerate(segs):
            phase = (codon_start - 1) if i == 0 else (3 - (done - (codon_start - 1)) % 3) % 3
            done += e - s + 1
            attrs = [("ID", f"{lt}.CDS"), ("Parent", mrna), ("Name", busco_id),
                     ("locus_tag", lt), ("product", product)] + common
            if i == 0:
                attrs += [("codon_start", codon_start),
                          ("translation_attached", "true" if "translation" in c.qualifiers else "false"),
                          ("Note", notes)]
            rows.append((s, 2, "\t".join([rec.id, SOURCE, "CDS", str(s), str(e), ".", strand,
                                          str(phase), _attrs(attrs)])))
    return rows, len(genes)


def build_gff3(records, fai_entries, scope: str) -> Tuple[str, Dict[str, int]]:
    order = [e[0] for e in fai_entries]
    lengths = {e[0]: e[1] for e in fai_entries}
    by_id = {r.id: r for r in records}
    if set(by_id) != set(order):
        missing = sorted(set(order) - set(by_id))[:5]
        extra = sorted(set(by_id) - set(order))[:5]
        raise ExportError(f"GenBank and FASTA contigs differ (only in FASTA: {missing}; "
                          f"only in GenBank: {extra})")
    lines = ["##gff-version 3"]
    lines += [f"##sequence-region {cid} 1 {lengths[cid]}" for cid in order]
    lines.append(f"# {scope}")
    counts = {}
    for cid in order:
        rows, n = gff_rows_for_record(by_id[cid])
        counts[cid] = n
        rows.sort(key=lambda t: (t[0], t[1], t[2]))
        lines += [t[2] for t in rows]
    return "\n".join(lines) + "\n", counts


def check_sequences(records, fasta: Path):
    """The FASTA must hold exactly the sequences of the GenBank records (no edits)."""
    gb = {r.id: r for r in records}
    n = 0
    for fr in SeqIO.parse(str(fasta), "fasta"):
        r = gb.get(fr.id)
        if r is None:
            raise ExportError(f"{fasta}: contig {fr.id} is not in the GenBank file")
        if str(r.seq).upper() != str(fr.seq).upper():
            raise ExportError(f"{fasta}: contig {fr.id} sequence differs from the GenBank record")
        n += 1
    if n != len(gb):
        raise ExportError(f"{fasta} has {n} contigs, the GenBank file {len(gb)}")


# ------------------------------------------------------------------ main ---
def find_tool(name: str, explicit: Optional[str] = None) -> Optional[str]:
    if not explicit:
        return shutil.which(name)
    found = shutil.which(explicit)
    if found:
        return found
    path = Path(explicit)
    if path.exists() and path.is_file():
        return str(path)
    return None


def require_tabix(explicit: Optional[str] = None) -> str:
    tabix = find_tool("tabix", explicit)
    if tabix:
        return tabix
    hint = f" at {explicit}" if explicit else " on PATH"
    raise ExportError(
        f"tabix not found{hint}. Install htslib/tabix or pass --tabix /path/to/tabix. "
        "Artemis needs the .tbi index for the whole-assembly GFF3 view."
    )


def export_sample(gbk: Path, fasta: Path, outdir: Path, stem: str,
                  tabix: Optional[str] = None) -> Dict:
    tabix = require_tabix(tabix)
    outdir.mkdir(parents=True, exist_ok=True)
    fasta_out = outdir / fasta.name
    if fasta_out.resolve() == fasta.resolve():
        raise ExportError("output directory must not be the FASTA's own directory")
    records = list(SeqIO.parse(str(gbk), "genbank"))
    if not records:
        raise ExportError(f"no records in {gbk}")
    scope_notes = [f.qualifiers.get("note", [""])[0] for f in records[0].features if f.type == "source"]
    scope = scope_notes[0] if scope_notes else "BuscoPainter: only BUSCO marker genes are annotated."

    entries = fasta_index(fasta)                      # fails early if faidx could not index it
    check_sequences(records, fasta)
    shutil.copyfile(fasta, fasta_out)
    if sha256(fasta_out) != sha256(fasta):
        raise ExportError(f"copy of {fasta} is not byte-identical")
    write_fai(entries, Path(str(fasta_out) + ".fai"))

    text, counts = build_gff3(records, entries, scope)
    gff_gz = outdir / f"{stem}.busco.gff3.gz"
    with bgzf.BgzfWriter(str(gff_gz), "wb") as w:
        w.write(text.encode("utf-8"))
    tbi = Path(str(gff_gz) + ".tbi")
    if tbi.exists():
        tbi.unlink()
    subprocess.run([tabix, "-p", "gff", str(gff_gz)], check=True)
    result = {"sample_stem": stem, "genbank": str(gbk), "fasta_source": str(fasta),
              "fasta_copy": str(fasta_out), "fasta_sha256": sha256(fasta_out),
              "gff3": str(gff_gz), "tbi": str(tbi) if tbi.exists() else None,
              "contigs": len(entries), "contigs_with_genes": sum(1 for v in counts.values() if v),
              "gene_models": sum(counts.values())}
    (outdir / "README.txt").write_text(
        "Open the whole assembly in Artemis\n\n"
        f"1. File > Open ...            {fasta_out.name}        -> click 'Use index'\n"
        f"2. File > Read An Entry ...   {gff_gz.name}   -> click 'Use index'\n"
        "3. Choose a contig in the drop-down list at the top of the window.\n\n"
        f"Terminal (if the art launcher is installed):  art {fasta_out.name} + {gff_gz.name}\n\n"
        + f"{scope}\nWritten by buscopainter-artemis {__version__} from {gbk.name}.\n")
    return result


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Write an indexed FASTA + indexed GFF3 pair so Artemis can browse every contig "
                    "of an assembly with its BuscoPainter annotations.")
    p.add_argument("--config", help="painter config (JSON); exports every sample")
    p.add_argument("--outdir", required=True,
                   help="with --config: the painter output directory (export goes to "
                        "<outdir>/<prefix>_<sample>_artemis/); with --gbk: the export directory")
    p.add_argument("--export-root", help="with --config: write the <prefix>_<sample>_artemis/ "
                   "folders here instead of in --outdir")
    p.add_argument("--prefix", help="painter prefix (default: config prefix, else Project)")
    p.add_argument("--samples", nargs="+")
    p.add_argument("--gbk", help="a single painter GenBank file (use with --fasta)")
    p.add_argument("--fasta", help="the assembly FASTA that GenBank file was made from")
    p.add_argument("--tabix", help="path to the tabix program (default: search PATH)")
    p.add_argument("--report", help="write a JSON summary")
    a = p.parse_args(argv)
    try:
        jobs = []
        if a.gbk or a.fasta:
            if not (a.gbk and a.fasta) or a.config:
                raise ExportError("use either --config, or --gbk together with --fasta")
            gbk = Path(a.gbk)
            jobs.append((gbk, Path(a.fasta), Path(a.outdir), gbk.stem))
        elif a.config:
            cfg = json.loads(Path(a.config).read_text())
            prefix = a.prefix or cfg.get("prefix", "Project")
            for s in cfg.get("samples", []):
                name = s.get("name") or Path(s["fasta"]).stem
                if a.samples and name not in a.samples:
                    continue
                pre = s.get("prefix") or prefix
                stem = f"{pre}_{name}"
                gbk = Path(a.outdir) / f"{stem}.gbk"
                if not gbk.exists():
                    raise ExportError(f"{gbk} not found; run buscopainter first (not --dry-run)")
                root = Path(a.export_root or a.outdir)
                jobs.append((gbk, Path(s["fasta"]), root / f"{stem}_artemis", stem))
            if not jobs:
                raise ExportError("no samples to export")
        else:
            raise ExportError("give --config or --gbk/--fasta")
        tabix = require_tabix(a.tabix)
        results = []
        for gbk, fasta, out, stem in jobs:
            r = export_sample(gbk, fasta, out, stem, tabix)
            print(f"[INFO] {stem}: {r['gene_models']} gene models on {r['contigs_with_genes']} of "
                  f"{r['contigs']} contigs -> {out}")
            results.append(r)
        if a.report:
            Path(a.report).write_text(json.dumps(results, indent=2))
        return 0
    except (ExportError, OSError, subprocess.CalledProcessError) as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
