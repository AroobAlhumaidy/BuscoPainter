#!/usr/bin/env python3
"""
Independent validation of BuscoPainter T7.1-fixed outputs.

Reads the ORIGINAL inputs (FASTA, BUSCO full_table.tsv, per-BUSCO GFFs, short_summary.json,
optionally BUSCO .faa) and the PRODUCED outputs (GenBank, skipped/issues/summary TSVs) and checks
them against each other. It deliberately does not import the pipeline; parsing, matching and
counting are re-implemented here in a different way (plain csv/text parsing, Biopython only for
FASTA/GenBank I/O and translation).

Usage:
  validate_busco_painter.py --config config.json --outdir outputs_t7_validation/run \
      [--prefix primary] [--samples TP2 ...] [--check-faa] [--check-miniprot] [--report out.json]
Exit code 0 only if every sample passes every check.
"""
import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

STATUS_SUBDIR = {"Complete": "single_copy_busco_sequences",
                 "Duplicated": "multi_copy_busco_sequences",
                 "Fragmented": "fragmented_busco_sequences"}


class Checks:
    def __init__(self, sample):
        self.sample = sample
        self.results = []  # (name, ok, detail)

    def add(self, name, ok, detail=""):
        self.results.append((name, bool(ok), str(detail)))

    @property
    def ok(self):
        return all(r[1] for r in self.results)


# ----------------------------------------------------------------- inputs ---
def find_run_dir(busco_dir: Path) -> Path:
    if (busco_dir / "full_table.tsv").exists():
        return busco_dir
    runs = [p for p in busco_dir.glob("run_*") if (p / "full_table.tsv").exists()]
    if len(runs) != 1:
        raise SystemExit(f"cannot find unique run_* in {busco_dir}")
    return runs[0]


def read_ft(path: Path):
    rows = []
    with open(path, newline="") as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            r = {"id": f[0], "status": f[1], "line": n}
            if f[1] != "Missing":
                a, b = int(f[3]), int(f[4])
                r.update(contig=f[2], lo=min(a, b), hi=max(a, b), strand=f[5], desc=f[9] if len(f) > 9 else "")
            rows.append(r)
    return rows


def read_gff_models(path: Path):
    """mRNA id -> dict(contig, start, end, strand, cds=[(start,end,phase)], stop=bool, attrs)."""
    models = {}
    kids = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        c = line.rstrip("\n").split("\t")
        attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
        if c[2] == "mRNA":
            models[attrs["ID"]] = {"contig": c[0], "start": int(c[3]), "end": int(c[4]), "strand": c[6],
                                   "cds": [], "stop": False, "attrs": attrs}
        else:
            kids.append((c, attrs))
    for c, attrs in kids:
        m = models[attrs["Parent"]]
        if c[2] == "CDS":
            m["cds"].append((int(c[3]), int(c[4]), int(c[7]), attrs))
        elif c[2] == "stop_codon":
            m["stop"] = True
    return models


def read_tsv(path: Path):
    if not path.exists():
        return None
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- checks ---
def validate_sample(s, outdir: Path, prefix: str, check_faa: bool, check_miniprot: bool = False):
    name = s["name"]
    ck = Checks(name)
    run = find_run_dir(Path(s["busco_dir"]))
    fasta = {r.id: str(r.seq) for r in SeqIO.parse(s["fasta"], "fasta")}
    ft = read_ft(run / "full_table.tsv")
    gbk_path = outdir / f"{prefix}_{name}.gbk"
    stats = {"sample": name}

    # --- GenBank readable --------------------------------------------------
    try:
        recs = list(SeqIO.parse(str(gbk_path), "genbank"))
        ck.add("genbank_parses", True, f"{len(recs)} records")
    except Exception as e:  # noqa: BLE001
        ck.add("genbank_parses", False, repr(e))
        return ck, stats
    gb = {r.id: r for r in recs}

    # --- contigs and sequences ---------------------------------------------
    ck.add("record_count_equals_fasta", len(recs) == len(fasta) == len(gb), f"gbk={len(recs)} fasta={len(fasta)}")
    ck.add("record_ids_equal_fasta_ids", set(gb) == set(fasta),
           f"missing={sorted(set(fasta) - set(gb))[:5]} extra={sorted(set(gb) - set(fasta))[:5]}")
    ck.add("record_order_equals_fasta", [r.id for r in recs] == list(fasta))
    bad_seq = [cid for cid in fasta if cid in gb and str(gb[cid].seq) != fasta[cid]]
    ck.add("sequences_identical", not bad_seq, f"{len(bad_seq)} differ {bad_seq[:5]}")
    stats["fasta_contigs"] = len(fasta)

    # --- features: bounds, strands, pairing --------------------------------
    genes, cdss = {}, {}
    oob, bad_strand, dup_tag, other_types = [], [], [], Counter()
    for r in recs:
        for f in r.features:
            if f.type == "source":
                if not (int(f.location.start) == 0 and int(f.location.end) == len(r)):
                    oob.append((r.id, "source"))
                continue
            if f.type not in ("gene", "CDS"):
                other_types[f.type] += 1
                continue
            parts = f.location.parts
            for p in parts:
                if not (0 <= int(p.start) < int(p.end) <= len(r)):
                    oob.append((r.id, f.qualifiers.get("locus_tag")))
            if {p.strand for p in parts} - {1, -1} or len({p.strand for p in parts}) != 1:
                bad_strand.append(f.qualifiers.get("locus_tag"))
            lt = f.qualifiers["locus_tag"][0]
            d = genes if f.type == "gene" else cdss
            if lt in d:
                dup_tag.append(lt)
            d[lt] = (r.id, f)
    ck.add("features_within_bounds", not oob, f"{len(oob)} {oob[:3]}")
    ck.add("feature_strands_valid", not bad_strand, bad_strand[:3])
    ck.add("locus_tags_unique_per_type", not dup_tag, dup_tag[:3])
    ck.add("only_source_gene_CDS_features", not other_types, dict(other_types))
    ck.add("gene_set_equals_CDS_set", set(genes) == set(cdss),
           f"genes_only={len(set(genes) - set(cdss))} cds_only={len(set(cdss) - set(genes))}")
    span_bad = []
    for lt, (cid, c) in cdss.items():
        if lt not in genes:
            continue
        g = genes[lt][1]
        if genes[lt][0] != cid or g.location.strand != c.location.strand or \
                int(g.location.start) > int(c.location.start) or int(g.location.end) < int(c.location.end):
            span_bad.append(lt)
    ck.add("CDS_inside_its_gene_same_strand", not span_bad, span_bad[:3])

    # --- trace each CDS to full_table + GFF --------------------------------
    present = [r for r in ft if r["status"] != "Missing"]
    missing_ids = {r["id"] for r in ft if r["status"] == "Missing"}
    rows_by_id = defaultdict(list)
    for r in present:
        rows_by_id[r["id"]].append(r)
    gff_cache = {}

    def models_for(bid, status):
        key = (bid, status)
        if key not in gff_cache:
            p = run / "busco_sequences" / STATUS_SUBDIR[status] / f"{bid}.gff"
            gff_cache[key] = read_gff_models(p) if p.exists() else {}
        return gff_cache[key]

    trace_fail, coord_fail, strand_fail, phase_fail, row_claims = [], [], [], [], defaultdict(list)
    missing_on_contig = []
    transl_bad, transl_n, multi_seg, trimmed = [], 0, 0, 0
    trimmed_lts = {}  # locus_tag -> GFF phase of the first exported segment
    for lt, (cid, c) in cdss.items():
        # locus_tag = <sample>_<buscoid>_<mRNA id>
        rest = lt[len(name) + 1:] if lt.startswith(name + "_") else None
        if rest is None or "_" not in rest:
            trace_fail.append((lt, "locus_tag format"))
            continue
        bid, mid = rest.rsplit("_", 1)
        if bid in missing_ids:
            missing_on_contig.append(lt)
        rows = rows_by_id.get(bid)
        if not rows:
            trace_fail.append((lt, "busco id not in full_table as present"))
            continue
        status = rows[0]["status"]
        m = models_for(bid, status).get(mid)
        if m is None:
            trace_fail.append((lt, f"mRNA {mid} not in {status} GFF"))
            continue
        if m["contig"] != cid:
            trace_fail.append((lt, "contig differs from GFF"))
            continue
        # full_table row on same contig/strand overlapping the model
        cand = [r for r in rows if r["contig"] == cid and r["lo"] <= m["end"] and m["start"] <= r["hi"]]
        if len(cand) != 1:
            # several rows may overlap; take the closest one and require uniqueness of the best
            cand.sort(key=lambda r: abs(r["lo"] - m["start"]) + abs(r["hi"] - m["end"]))
            if not cand or (len(cand) > 1 and abs(cand[0]["lo"] - m["start"]) + abs(cand[0]["hi"] - m["end"]) ==
                            abs(cand[1]["lo"] - m["start"]) + abs(cand[1]["hi"] - m["end"])):
                trace_fail.append((lt, f"{len(cand)} full_table rows match"))
                continue
        row = cand[0]
        row_claims[row["line"]].append(lt)
        strand_num = 1 if m["strand"] == "+" else -1
        if c.location.strand != strand_num or row["strand"] != m["strand"]:
            strand_fail.append(lt)
        # Expected CDS rows = BUSCO's prediction. If the full_table span is a run of whole GFF CDS
        # rows (others wholly outside), BUSCO dropped the outer rows and so must the GenBank;
        # otherwise every GFF CDS row is expected.
        rows_exp = m["cds"]
        if abs(row["lo"] - m["start"]) > 3 or abs(row["hi"] - m["end"]) > 3:
            inside = [x for x in m["cds"] if x[0] >= row["lo"] - 3 and x[1] <= row["hi"] + 3]
            outside = [x for x in m["cds"] if x[1] < row["lo"] or x[0] > row["hi"]]
            if inside and len(inside) + len(outside) == len(m["cds"]) and \
                    abs(min(x[0] for x in inside) - row["lo"]) <= 3 and abs(max(x[1] for x in inside) - row["hi"]) <= 3:
                rows_exp = inside
                trimmed += 1
                trimmed_lts[lt] = sorted(inside, key=lambda x: x[0], reverse=(m["strand"] == "-"))[0][2]
        # gene == span of the expected CDS rows (= mRNA span when nothing is dropped)
        g = genes[lt][1]
        g_lo = m["start"] if rows_exp is m["cds"] else min(x[0] for x in rows_exp)
        g_hi = m["end"] if rows_exp is m["cds"] else max(x[1] for x in rows_exp)
        if (int(g.location.start), int(g.location.end)) != (g_lo - 1, g_hi):
            coord_fail.append((lt, "gene span"))
        want = sorted(((a - 1, b) for a, b, _, _ in rows_exp), reverse=(strand_num == -1))
        got = [(int(p.start), int(p.end)) for p in c.location.parts]
        if got != want:
            coord_fail.append((lt, f"CDS parts {got[:3]} != GFF {want[:3]}"))
        first_phase = sorted(rows_exp, key=lambda x: x[0], reverse=(strand_num == -1))[0][2]
        if c.qualifiers.get("codon_start") != [str(first_phase + 1)]:
            phase_fail.append(lt)
        if len(rows_exp) > 1:
            multi_seg += 1
        # translation (if attached) must be the exact clean translation of the CDS
        if "translation" in c.qualifiers:
            transl_n += 1
            nt = c.location.extract(gb[cid].seq)[first_phase:]
            p = str(nt.translate(table=int(c.qualifiers.get("transl_table", ["1"])[0])))
            p = p[:-1] if p.endswith("*") else p
            if len(nt) % 3 or "*" in p or p != c.qualifiers["translation"][0]:
                transl_bad.append(lt)
    ck.add("every_feature_traced_to_busco_and_gff", not trace_fail, f"{len(trace_fail)} {trace_fail[:3]}")
    ck.add("gene_and_CDS_coordinates_equal_GFF", not coord_fail, f"{len(coord_fail)} {coord_fail[:3]}")
    ck.add("strand_matches_GFF_and_full_table", not strand_fail, strand_fail[:3])
    ck.add("codon_start_equals_first_segment_phase", not phase_fail, phase_fail[:3])
    ck.add("attached_translations_valid", not transl_bad, f"{len(transl_bad)} of {transl_n} {transl_bad[:3]}")
    ck.add("missing_buscos_not_annotated", not missing_on_contig, missing_on_contig[:3])
    double = {k: v for k, v in row_claims.items() if len(v) > 1}
    ck.add("no_full_table_row_used_twice", not double, list(double.items())[:3])

    # --- reconciliation with full_table and skipped report -----------------
    skipped = read_tsv(outdir / f"{prefix}_{name}.skipped.tsv") or []
    sk_keys = Counter((r["busco_id"], r["sequence"]) for r in skipped)
    unexplained = []
    claimed_lines = set(row_claims)
    for r in present:
        if r["line"] in claimed_lines:
            continue
        k = (r["id"], r["contig"])
        if sk_keys[k] > 0:
            sk_keys[k] -= 1
        else:
            unexplained.append((r["id"], r["contig"], r["lo"], r["hi"]))
    ck.add("every_full_table_hit_exported_or_reported", not unexplained, f"{len(unexplained)} {unexplained[:3]}")
    phantom = [k for k, v in sk_keys.items() if v > 0]
    ck.add("skipped_rows_refer_to_real_hits", not phantom, phantom[:3])

    # duplicated: distinct occurrences retained
    dup_ids = {r["id"] for r in present if r["status"] == "Duplicated"}
    dup_loss = []
    for bid in dup_ids:
        rows = rows_by_id[bid]
        exported = [cdss[lt] for lts in (row_claims.get(r["line"], []) for r in rows) for lt in lts]
        locs = {(cid, int(f.location.start), int(f.location.end)) for cid, f in exported}
        n_skip = sum(1 for r in skipped if r["busco_id"] == bid)
        if len(locs) != len(exported) or len(exported) + n_skip != len(rows):
            dup_loss.append(bid)
    ck.add("duplicated_copies_distinct_and_complete", not dup_loss, f"{len(dup_loss)} {dup_loss[:3]}")

    # counts: ortholog-level vs hit-level, against short_summary.json and pipeline summary
    by_status_ids = defaultdict(set)
    by_status_rows = Counter()
    for r in ft:
        by_status_ids[r["status"]].add(r["id"])
        by_status_rows[r["status"]] += 1
    S, D, F, M = (len(by_status_ids[k]) for k in ("Complete", "Duplicated", "Fragmented", "Missing"))
    stats.update(S_orthologs=S, D_orthologs=D, C_orthologs=S + D, F_orthologs=F, M_orthologs=M,
                 S_rows=by_status_rows["Complete"], D_rows=by_status_rows["Duplicated"],
                 F_rows=by_status_rows["Fragmented"], hit_rows=len(present),
                 genes=len(genes), CDS=len(cdss),
                 cds_segments=sum(len(c.location.parts) for _, c in cdss.values()),
                 multi_segment=multi_seg, span_trimmed=trimmed, translations=transl_n, rejected=len(skipped))
    ssj = run / "short_summary.json"
    if ssj.exists():
        res = json.load(open(ssj))["results"]
        exp = (res["Single copy BUSCOs"], res["Multi copy BUSCOs"], res["Fragmented BUSCOs"],
               res["Missing BUSCOs"], res["n_markers"])
        got = (S, D, F, M, len({r["id"] for r in ft}))
        ck.add("full_table_ortholog_counts_equal_short_summary", tuple(map(int, exp)) == got, f"json={exp} ft={got}")
    summ = read_tsv(outdir / f"{prefix}.busco_summary.tsv")
    srow = next((r for r in summ or [] if r["sample"] == name), None)
    if srow:
        pairs = {"busco_S_orthologs": S, "busco_D_orthologs": D, "busco_F_orthologs": F, "busco_M_orthologs": M,
                 "busco_C_orthologs": S + D, "busco_S_hit_rows": by_status_rows["Complete"],
                 "busco_D_hit_rows": by_status_rows["Duplicated"], "busco_F_hit_rows": by_status_rows["Fragmented"],
                 "exported_genes": len(genes), "exported_CDS": len(cdss), "rejected_hits": len(skipped),
                 "fasta_contigs": len(fasta)}
        diff = {k: (srow.get(k), v) for k, v in pairs.items() if str(srow.get(k)) != str(v)}
        ck.add("pipeline_summary_matches_independent_counts", not diff, diff)
    else:
        ck.add("pipeline_summary_matches_independent_counts", False, "sample row missing")
    ck.add("exported_plus_rejected_equals_hit_rows", len(cdss) + len(skipped) == len(present),
           f"{len(cdss)}+{len(skipped)} vs {len(present)}")

    # --- optional cross-check of attached translations with BUSCO .faa -----
    if check_faa:
        cat = Counter()
        examples = []
        faa_cache = {}
        for lt, (cid, c) in cdss.items():
            if "translation" not in c.qualifiers:
                continue
            bid, mid = lt[len(name) + 1:].rsplit("_", 1)
            rows = rows_by_id[bid]
            status = rows[0]["status"]
            fp = run / "busco_sequences" / STATUS_SUBDIR[status] / f"{bid}.faa"
            if fp not in faa_cache:
                faa_cache[fp] = {r.id.split("|", 1)[1]: str(r.seq) for r in SeqIO.parse(str(fp), "fasta")} if fp.exists() else {}
            m = models_for(bid, status)[mid]
            row = next(r for r in rows if row_claims.get(r["line"]) == [lt])
            bus = faa_cache[fp].get(f"{cid}:{row['lo']}-{row['hi']}|{row['strand']}")
            ours = c.qualifiers["translation"][0]
            if bus is None:
                cat["faa_record_not_found"] += 1
                continue
            b = bus.rstrip("*")
            if b == ours:
                cat["identical"] += 1
            elif not m["stop"] and b[:-1] == ours:
                cat["busco_faa_has_one_extra_C_terminal_residue_no_stop_model"] += 1
            elif trimmed_lts.get(lt) and ours in (b[1:], b[1:-1]):
                # BUSCO keeps the residue whose codon spans the dropped/kept segment junction;
                # its first bases lie outside the exported CDS (codon_start > 1)
                cat["busco_faa_has_one_extra_N_terminal_junction_residue_trimmed_model"] += 1
            else:
                cat["different"] += 1
                examples.append(lt)
        stats["faa_crosscheck"] = dict(cat)
        ck.add("translations_agree_with_busco_faa", cat["different"] == 0 and cat["faa_record_not_found"] == 0,
               f"{dict(cat)} {examples[:3]}")

    # --- optional cross-check with miniprot's own translation (##STA) ------
    if check_miniprot:
        logs = sorted((Path(s["busco_dir"]) / "logs").glob("miniprot_align_*_out.log"))
        if len(logs) != 1:
            ck.add("translations_agree_with_miniprot_STA", False, f"{len(logs)} miniprot logs found")
        else:
            sta, prev = {}, None
            with open(logs[0]) as fh:
                for line in fh:
                    if line.startswith("##STA\t"):
                        prev = line.rstrip("\n").split("\t", 1)[1]
                    elif "\tmRNA\t" in line and prev is not None:
                        mid = line.split("ID=", 1)[1].split(";", 1)[0]
                        sta[mid] = prev
                        prev = None
            cat, examples = Counter(), []
            for lt, (cid, c) in cdss.items():
                if "translation" not in c.qualifiers:
                    continue
                bid, mid = lt[len(name) + 1:].rsplit("_", 1)
                p = sta.get(mid)
                ours = c.qualifiers["translation"][0]
                has_stop = models_for(bid, rows_by_id[bid][0]["status"])[mid]["stop"]
                if p is None:
                    cat["STA_not_found"] += 1
                elif p.rstrip("*") == ours:
                    cat["identical"] += 1
                elif not has_stop and p[:-1] == ours:
                    # miniprot ##STA appends the codon after the aligned region when the
                    # alignment has no stop codon (GFF CDS = PAF span = CIGAR length)
                    cat["STA_has_one_extra_C_terminal_residue_no_stop_model"] += 1
                elif lt in trimmed_lts and ours in p:
                    # ##STA is the untrimmed model; the BUSCO-trimmed CDS is a part of it
                    cat["trimmed_model_translation_inside_STA"] += 1
                else:
                    cat["different"] += 1
                    examples.append(lt)
            stats["miniprot_STA_crosscheck"] = dict(cat)
            ck.add("translations_agree_with_miniprot_STA", cat["different"] == 0 and cat["STA_not_found"] == 0,
                   f"{dict(cat)} {examples[:3]}")

    # --- issue table summary (for the report) ------------------------------
    iss = read_tsv(outdir / f"{prefix}_{name}.annotation_issues.tsv") or []
    stats["issues"] = dict(Counter(f"{r['severity']}:{r['code']}" for r in iss))
    return ck, stats


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix")
    ap.add_argument("--samples", nargs="+")
    ap.add_argument("--check-faa", action="store_true", help="compare translations with BUSCO .faa")
    ap.add_argument("--check-miniprot", action="store_true", help="compare translations with miniprot ##STA in logs/")
    ap.add_argument("--report", help="write JSON report here")
    a = ap.parse_args(argv)
    cfg = json.load(open(a.config))
    prefix = a.prefix or cfg.get("prefix", "Project")
    outdir = Path(a.outdir)
    samples = [s for s in cfg["samples"] if not a.samples or s["name"] in a.samples]
    report, all_ok = [], True
    for s in samples:
        ck, stats = validate_sample(s, outdir, prefix, a.check_faa, a.check_miniprot)
        all_ok &= ck.ok
        print(f"== {s['name']}: {'PASS' if ck.ok else 'FAIL'}")
        for n, ok, d in ck.results:
            print(f"   [{'ok' if ok else 'FAIL'}] {n}" + (f"  ({d})" if (not ok or d) and len(d) < 300 else ""))
        print("   stats:", json.dumps(stats))
        report.append({"sample": s["name"], "pass": ck.ok, "stats": stats,
                       "checks": [{"check": n, "ok": ok, "detail": d} for n, ok, d in ck.results]})
    if a.report:
        Path(a.report).write_text(json.dumps(report, indent=2))
    print("OVERALL:", "PASS" if all_ok else "FAIL")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
