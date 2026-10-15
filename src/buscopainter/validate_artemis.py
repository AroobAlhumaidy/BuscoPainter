#!/usr/bin/env python3
"""
Independent check of a ``buscopainter-artemis`` export (indexed FASTA + indexed GFF3).

It does not import the exporter. For each sample it re-reads the source FASTA, the painter
GenBank and the painter ``features.tsv`` and checks the export directory:

  fasta_copy_identical   the exported FASTA is byte-identical to the source FASTA
  fai_valid              the .fai matches an index recomputed here (and samtools, if found)
  bgzf_valid             the GFF3 is BGZF (every block has the BC header, EOF block present)
  sequence_regions       one ##sequence-region per FASTA contig, with the right length
  gff3_structure         9 columns, valid strand/phase, unique IDs, Parents exist, sorted
  models_vs_genbank      every locus_tag: same contig, gene span, strand and CDS segments
  models_vs_features_tsv the same against features.tsv (busco id, status, codon_start too)
  phases_consistent      each CDS row's phase follows from codon_start and segment lengths
  counts_per_contig      gene models per contig identical in GFF3, GenBank and features.tsv
  tabix_index            .tbi present; with `tabix`, each contig query returns its rows

Usage:
  buscopainter-validate-artemis --config config.json --outdir results [--report out.json]
Exit 0 only if every check passes for every sample.
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO

BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def recompute_fai(fasta: Path):
    out, name, length, off, lb, lw = [], None, 0, 0, None, None
    pos = 0
    with open(fasta, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    out.append(f"{name}\t{length}\t{off}\t{lb}\t{lw}")
                name, length, lb, lw = line[1:].split()[0].decode(), 0, None, None
                off = pos + len(line)
            else:
                s = line.rstrip(b"\r\n")
                if s and lb is None:
                    lb, lw = len(s), len(line)
                length += len(s)
            pos += len(line)
    if name is not None:
        out.append(f"{name}\t{length}\t{off}\t{lb}\t{lw}")
    return "\n".join(out) + "\n"


def bgzf_ok(path: Path):
    data = path.read_bytes()
    i, blocks = 0, 0
    while i < len(data):
        if data[i:i + 4] != b"\x1f\x8b\x08\x04" or data[i + 12:i + 14] != b"BC":
            return False, f"block {blocks} at byte {i} is not a BGZF block"
        bsize = int.from_bytes(data[i + 16:i + 18], "little") + 1
        i += bsize
        blocks += 1
    if not data.endswith(BGZF_EOF):
        return False, "BGZF EOF marker missing"
    return True, f"{blocks} blocks"


def parse_attrs(col9):
    out = {}
    for kv in col9.split(";"):
        if kv:
            k, v = kv.split("=", 1)
            out[k] = [unquote(x) for x in v.split(",")]
    return out


def read_gff3(path: Path):
    regions, rows, problems = {}, [], []
    with gzip.open(path, "rt") as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if line.startswith("##sequence-region"):
                _, sid, s, e = line.split()
                regions[sid] = (int(s), int(e))
                continue
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) != 9:
                problems.append(f"line {n}: {len(f)} columns")
                continue
            rows.append({"n": n, "seqid": f[0], "type": f[2], "start": int(f[3]), "end": int(f[4]),
                         "strand": f[6], "phase": f[7], "a": parse_attrs(f[8])})
    return regions, rows, problems


def genbank_models(gbk: Path):
    models, lengths = {}, {}
    for rec in SeqIO.parse(str(gbk), "genbank"):
        lengths[rec.id] = len(rec)
        for f in rec.features:
            if f.type in ("gene", "CDS"):
                lt = f.qualifiers["locus_tag"][0]
                m = models.setdefault(lt, {"contig": rec.id})
                m["strand"] = "+" if f.location.strand == 1 else "-"
                segs = sorted((int(p.start) + 1, int(p.end)) for p in f.location.parts)
                m[f.type] = segs
                if f.type == "CDS":
                    m["codon_start"] = int(f.qualifiers.get("codon_start", ["1"])[0])
    return models, lengths


def segs_from_loc_string(loc: str):
    return sorted((int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", loc))


def check_sample(name, fasta, gbk, features_tsv, exp_dir, stem, tabix):
    res = []

    def add(check, ok, detail=""):
        res.append({"check": check, "ok": bool(ok), "detail": str(detail)[:500]})

    fa_copy = exp_dir / fasta.name
    gff = exp_dir / f"{stem}.busco.gff3.gz"
    for p in (fa_copy, Path(str(fa_copy) + ".fai"), gff):
        if not p.exists():
            add("files_present", False, f"missing {p}")
            return res
    add("files_present", True)
    add("fasta_copy_identical", _sha(fa_copy) == _sha(fasta), fa_copy)

    fai_txt = Path(str(fa_copy) + ".fai").read_text()
    ok = fai_txt == recompute_fai(fasta)
    detail = "recomputed"
    samtools = shutil.which("samtools")
    if samtools and ok:
        tmp = exp_dir / ".validate_tmp.fasta"
        shutil.copyfile(fasta, tmp)
        try:
            subprocess.run([samtools, "faidx", str(tmp)], check=True)
            ok = Path(str(tmp) + ".fai").read_text() == fai_txt
            detail = "recomputed + samtools faidx"
        finally:
            for t in (tmp, Path(str(tmp) + ".fai")):
                t.unlink(missing_ok=True)
    add("fai_valid", ok, detail)

    ok, detail = bgzf_ok(gff)
    add("bgzf_valid", ok, detail)

    fai = [l.split("\t") for l in fai_txt.splitlines()]
    fai_len = {r[0]: int(r[1]) for r in fai}
    regions, rows, problems = read_gff3(gff)
    add("sequence_regions", regions == {k: (1, v) for k, v in fai_len.items()},
        f"{len(regions)} regions, {len(fai_len)} contigs")

    # structure
    ids = Counter()
    by_id = defaultdict(list)
    order = [r[0] for r in fai]
    last = (-1, 0)
    for r in rows:
        if r["seqid"] not in fai_len:
            problems.append(f"line {r['n']}: unknown seqid {r['seqid']}")
            continue
        if not (1 <= r["start"] <= r["end"] <= fai_len[r["seqid"]]):
            problems.append(f"line {r['n']}: out of range")
        if r["strand"] not in "+-" or r["phase"] not in (("0", "1", "2") if r["type"] == "CDS" else (".",)):
            problems.append(f"line {r['n']}: strand/phase {r['strand']}/{r['phase']}")
        key = (order.index(r["seqid"]), r["start"])
        if key < last:
            problems.append(f"line {r['n']}: not sorted")
        last = key
        rid = r["a"].get("ID", [None])[0]
        by_id[rid].append(r)
        if r["type"] != "CDS":
            ids[rid] += 1
    problems += [f"duplicate ID {k}" for k, v in ids.items() if v > 1]
    for r in rows:
        for p in r["a"].get("Parent", []):
            if p not in by_id:
                problems.append(f"line {r['n']}: Parent {p} not found")
    add("gff3_structure", not problems, "; ".join(problems[:5]) or f"{len(rows)} rows")

    # models from the GFF3
    gm = {}
    for r in rows:
        lt = r["a"]["locus_tag"][0]
        m = gm.setdefault(lt, {"CDS": [], "gene": None, "mRNA": None})
        if r["type"] == "CDS":
            m["CDS"].append(r)
        else:
            m[r["type"]] = r
    gb, gb_len = genbank_models(gbk)
    bad = []
    for lt, w in gb.items():
        m = gm.get(lt)
        if not m or not m["gene"] or not m["mRNA"] or not m["CDS"]:
            bad.append(f"{lt}: missing rows")
            continue
        g = m["gene"]
        if g["seqid"] != w["contig"] or [(g["start"], g["end"])] != w["gene"] or g["strand"] != w["strand"]:
            bad.append(f"{lt}: gene differs")
        if (m["mRNA"]["start"], m["mRNA"]["end"]) != (g["start"], g["end"]):
            bad.append(f"{lt}: mRNA span differs from gene")
        segs = sorted((c["start"], c["end"]) for c in m["CDS"])
        if segs != w["CDS"] or {c["strand"] for c in m["CDS"]} != {w["strand"]}:
            bad.append(f"{lt}: CDS {segs} vs GenBank {w['CDS']}")
        if len({c['a']['ID'][0] for c in m["CDS"]}) != 1:
            bad.append(f"{lt}: CDS rows have different IDs")
    extra = set(gm) - set(gb)
    add("models_vs_genbank", not bad and not extra and gb_len == fai_len,
        "; ".join(bad[:5] + [f"extra {sorted(extra)[:3]}"] * bool(extra)
                  + ["GenBank/FASTA contigs differ"] * (gb_len != fai_len)) or f"{len(gb)} models")

    # features.tsv: a third, tabular source written by the painter
    bad = []
    with open(features_tsv, newline="") as fh:
        ft = {r["locus_tag"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for lt, r in ft.items():
        m = gm.get(lt)
        if not m or not m["CDS"]:
            bad.append(f"{lt}: not in GFF3")
            continue
        first = min(m["CDS"], key=lambda c: c["start"] if r["strand"] == "+" else -c["end"])
        a = first["a"]
        if (a.get("busco_id", [""])[0] != r["busco_id"] or a.get("busco_status", [""])[0] != r["status"]
                or m["gene"]["seqid"] != r["contig"] or m["gene"]["strand"] != r["strand"]
                or (m["gene"]["start"], m["gene"]["end"]) != (int(r["gene_start"]), int(r["gene_end"]))
                or sorted((c["start"], c["end"]) for c in m["CDS"]) != segs_from_loc_string(r["cds_location"])
                or a.get("codon_start", [""])[0] != r["codon_start"]):
            bad.append(f"{lt}: differs from features.tsv")
    add("models_vs_features_tsv", not bad and set(ft) == set(gm),
        "; ".join(bad[:5]) or f"{len(ft)} rows")

    # phases
    bad = []
    for lt, m in gm.items():
        cds = sorted(m["CDS"], key=lambda c: c["start"], reverse=m["CDS"][0]["strand"] == "-")
        cs = int(cds[0]["a"]["codon_start"][0])
        done = 0
        for i, c in enumerate(cds):
            want = cs - 1 if i == 0 else (-(done - (cs - 1))) % 3
            if int(c["phase"]) != want:
                bad.append(f"{lt} segment {i}: phase {c['phase']} expected {want}")
            done += c["end"] - c["start"] + 1
    add("phases_consistent", not bad, "; ".join(bad[:5]))

    gff_counts = Counter(m["gene"]["seqid"] for m in gm.values() if m["gene"])
    gb_counts = Counter(w["contig"] for w in gb.values())
    ft_counts = Counter(r["contig"] for r in ft.values())
    add("counts_per_contig", gff_counts == gb_counts == ft_counts,
        f"{sum(gff_counts.values())} models on {len(gff_counts)} contigs")

    tbi = Path(str(gff) + ".tbi")
    if not tbi.exists():
        add("tabix_index", False, "no .tbi: Artemis cannot use this GFF3 (run tabix -p gff)")
    elif tabix:
        listed = subprocess.run([tabix, "-l", str(gff)], capture_output=True, text=True, check=True).stdout.split()
        rows_per = Counter(r["seqid"] for r in rows)
        bad = [] if sorted(listed) == sorted(rows_per) else ["tabix -l contigs differ"]
        for sid in listed:
            out = subprocess.run([tabix, str(gff), sid], capture_output=True, text=True, check=True).stdout
            if len(out.splitlines()) != rows_per[sid]:
                bad.append(f"{sid}: tabix {len(out.splitlines())} rows, file {rows_per[sid]}")
        add("tabix_index", not bad, "; ".join(bad[:5]) or f"{len(listed)} contigs queried")
    else:
        add("tabix_index", True, ".tbi present; tabix not found, queries not checked")
    return res


def main(argv=None):
    p = argparse.ArgumentParser(description="Validate buscopainter-artemis exports.")
    p.add_argument("--config", required=True)
    p.add_argument("--outdir", required=True, help="painter output directory")
    p.add_argument("--export-root", help="where the <prefix>_<sample>_artemis/ folders are "
                   "(default: --outdir)")
    p.add_argument("--prefix")
    p.add_argument("--samples", nargs="+")
    p.add_argument("--tabix", help="path to tabix (default: search PATH)")
    p.add_argument("--report")
    a = p.parse_args(argv)
    cfg = json.loads(Path(a.config).read_text())
    prefix = a.prefix or cfg.get("prefix", "Project")
    tabix = a.tabix or shutil.which("tabix")
    out, all_ok = {}, True
    for s in cfg.get("samples", []):
        name = s.get("name") or Path(s["fasta"]).stem
        if a.samples and name not in a.samples:
            continue
        stem = f"{s.get('prefix') or prefix}_{name}"
        od = Path(a.outdir)
        res = check_sample(name, Path(s["fasta"]), od / f"{stem}.gbk", od / f"{stem}.features.tsv",
                           Path(a.export_root or od) / f"{stem}_artemis", stem, tabix)
        out[name] = res
        ok = all(r["ok"] for r in res)
        all_ok &= ok
        print(f"== {name}: {'PASS' if ok else 'FAIL'}")
        for r in res:
            print(f"  [{'ok' if r['ok'] else 'FAIL'}] {r['check']}: {r['detail']}")
    if not out:
        print("no samples checked")
        all_ok = False
    print(f"OVERALL: {'PASS' if all_ok else 'FAIL'}")
    if a.report:
        Path(a.report).write_text(json.dumps(out, indent=2))
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
