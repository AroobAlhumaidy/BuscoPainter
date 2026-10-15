#!/usr/bin/env python3
"""Compare what the real Artemis GUI loaded (ArtemisGuiCheck output) with the T7.1-fixed.3 GenBank.

For every contig in the GenBank: did Artemis offer it in the drop-down, with the right length,
and does it hold the same gene models (locus_tag), each with identical gene span, strand and CDS
segments (joins)? Also checks what the feature view drew while scrolling (display.tsv).

Usage: compare_artemis_vs_genbank.py <genbank> <ArtemisGuiCheck outdir> [--fasta-entry NAME]
Prints a summary and writes <outdir>/comparison.json; exit 0 only if everything matches.
"""
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


def segs_from_artemis(loc: str):
    # Artemis writes a minus-strand join either as complement(join(a,b)) or as
    # join(complement(b),complement(a)); both mean every part is on the minus strand.
    pairs = [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", loc)]
    n_comp = loc.count("complement(")
    if loc.startswith("complement("):
        strand = "-" if n_comp == 1 else "?"
    else:
        strand = "-" if n_comp and n_comp == len(pairs) else ("+" if n_comp == 0 else "?")
    if not pairs:  # single-base location
        pairs = [(int(x), int(x)) for x in re.findall(r"\d+", loc)[:1]]
    return strand, tuple(sorted(pairs))


def genbank_models(gbk: Path):
    contigs, models = {}, {}
    for rec in SeqIO.parse(str(gbk), "genbank"):
        contigs[rec.id] = len(rec)
        for f in rec.features:
            if f.type not in ("gene", "CDS"):
                continue
            lt = f.qualifiers["locus_tag"][0]
            m = models.setdefault(lt, {"contig": rec.id})
            strand = "+" if f.location.strand == 1 else "-"
            segs = tuple(sorted((int(p.start) + 1, int(p.end)) for p in f.location.parts))
            m[f.type] = (strand, segs)
            if f.type == "CDS":
                m["busco_id"] = re.match(r"BUSCO (\S+) ", f.qualifiers["note"][0]).group(1)
    return contigs, models


def main(gbk, outdir):
    out = Path(outdir)
    contigs, gb = genbank_models(Path(gbk))
    rows = list(csv.DictReader(open(out / "contigs.tsv"), delimiter="\t"))
    feats = list(csv.DictReader(open(out / "features.tsv"), delimiter="\t"))
    disp = list(csv.DictReader(open(out / "display.tsv"), delimiter="\t"))
    res = {"genbank_contigs": len(contigs), "dropdown_contigs": len(rows)}
    labels = [r["label"] for r in rows]
    res["dropdown_matches_genbank_order"] = labels == list(contigs)
    res["length_mismatches"] = [r["label"] for r in rows
                                if contigs.get(r["label"]) != int(r["artemis_length"])]

    # Artemis features grouped by locus_tag and key
    art = defaultdict(lambda: defaultdict(list))
    for f in feats:
        if f["locus_tag"]:
            art[f["locus_tag"]][f["key"]].append((f["label"], segs_from_artemis(f["location"])))
    problems = defaultdict(list)
    for lt, m in gb.items():
        a = art.get(lt)
        if not a:
            problems["missing_in_artemis"].append(lt)
            continue
        for key in ("gene", "CDS"):
            got = a.get(key, [])
            want = m[key]
            if len(got) != 1:
                problems[f"{key}_count_not_1"].append(f"{lt}:{len(got)}")
                continue
            contig, (strand, segs) = got[0]
            if contig != m["contig"]:
                problems[f"{key}_wrong_contig"].append(lt)
            if strand != want[0]:
                problems[f"{key}_strand"].append(lt)
            if segs != want[1]:
                problems[f"{key}_coords_or_joins"].append(f"{lt}: artemis {segs} genbank {want[1]}")
        if len(a.get("mRNA", [])) != 1:
            problems["mRNA_count_not_1"].append(lt)
    res["extra_in_artemis"] = sorted(set(art) - set(gb))[:20]
    res["genbank_models"] = len(gb)
    res["artemis_models"] = len(set(art) & set(gb))
    res["multi_segment_models"] = sum(1 for m in gb.values() if len(m["CDS"][1]) > 1)
    res["problems"] = {k: {"n": len(v), "examples": v[:5]} for k, v in problems.items()}

    # what the feature view drew while scrolling (CDS keyed by locus_tag)
    drawn = {d["locus_tag"] for d in disp if d["key"] == "CDS" and d["locus_tag"]}
    res["drawn_CDS_models"] = len(drawn & set(gb))
    res["not_drawn"] = sorted(set(gb) - drawn)[:20]
    per_contig = defaultdict(int)
    for m in gb.values():
        per_contig[m["contig"]] += 1
    res["contigs_with_models"] = len(per_contig)
    res["PASS"] = (res["dropdown_matches_genbank_order"] and not res["length_mismatches"]
                   and not problems and not res["extra_in_artemis"]
                   and res["drawn_CDS_models"] == len(gb))
    (out / "comparison.json").write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k != "problems"}, indent=1))
    for k, v in res["problems"].items():
        print(f"PROBLEM {k}: {v['n']}  e.g. {v['examples'][:2]}")
    print("PASS" if res["PASS"] else "FAIL")
    return 0 if res["PASS"] else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
