"""Edge-case tests for buscopainter.painter (reference: buscoPainter.T7.1.fixed.py) using a
small synthetic BUSCO run.

Every expected value below is derived from the synthetic sequences built here, not from the
pipeline's own helpers. Run with:  python -m pytest -q
"""
import json
import random
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from buscopainter import painter as bp

SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")]


def orf(n_codons, seed, internal_stop_at=None):
    rnd = random.Random(seed)
    cod = [rnd.choice(SENSE) for _ in range(n_codons)]
    if internal_stop_at is not None:
        cod[internal_stop_at] = "TGA"
    return "ATG" + "".join(cod) + "TAA"


def rc(s):
    return str(Seq(s).reverse_complement())


class Dataset:
    """Builds a FASTA + BUSCO v6-like run dir with full_table.tsv and per-BUSCO GFFs."""

    def __init__(self, root: Path):
        self.root = root
        rnd = random.Random(1)
        self.contigs = {"chrA": [rnd.choice("ACGT") for _ in range(3000)],
                        "chrB": [rnd.choice("ACGT") for _ in range(600)]}
        self.ft = []
        self.gff = {}  # (status_dir, busco_id) -> lines
        self.n = 0

    def place(self, contig, start1, seq):
        self.contigs[contig][start1 - 1:start1 - 1 + len(seq)] = list(seq)

    def model(self, busco_id, status, contig, strand, cds, stop=None, attrs="", sdir=None):
        """cds: [(start, end, phase)] 1-based inclusive. Appends the GFF mRNA and its rows."""
        self.n += 1
        mid = f"MP{self.n:04d}"
        lo, hi = min(c[0] for c in cds), max(c[1] for c in cds)
        a = f"ID={mid};Rank=1;Identity=0.9;Target={busco_id}_x 1 99" + (f";{attrs}" if attrs else "")
        lines = [f"{contig}\tminiprot\tmRNA\t{lo}\t{hi}\t100\t{strand}\t.\t{a}"]
        for s, e, ph in cds:
            lines.append(f"{contig}\tminiprot\tCDS\t{s}\t{e}\t100\t{strand}\t{ph}\tParent={mid};Rank=1;Target=x 1 9")
        if stop:
            lines.append(f"{contig}\tminiprot\tstop_codon\t{stop[0]}\t{stop[1]}\t0\t{strand}\t0\tParent={mid};Rank=1")
        sdir = sdir or bp.STATUS_DIRS[status]
        self.gff.setdefault((sdir, busco_id), []).extend(lines)
        return mid, lo, hi

    def row(self, busco_id, status, contig="", lo=None, hi=None, strand="", desc="test protein"):
        if status == "Missing":
            self.ft.append(f"{busco_id}\tMissing")
            return
        # BUSCO convention: start/end written in strand order
        a, b = (lo, hi) if strand == "+" else (hi, lo)
        self.ft.append(f"{busco_id}\t{status}\t{contig}\t{a}\t{b}\t{strand}\t100.0\t50\turl\t{desc}")

    def write(self, extra_cfg=None):
        fa = self.root / "asm.fasta"
        with open(fa, "w") as fh:
            for cid, s in self.contigs.items():
                fh.write(f">{cid}\n{''.join(s)}\n")
        run = self.root / "busco" / "run_test_odb12"
        for (sdir, bid), lines in self.gff.items():
            d = run / "busco_sequences" / sdir
            d.mkdir(parents=True, exist_ok=True)
            (d / f"{bid}.gff").write_text("\n".join(lines) + "\n")
        for sdir in bp.STATUS_DIRS.values():
            (run / "busco_sequences" / sdir).mkdir(parents=True, exist_ok=True)
        hdr = "# BUSCO version is: 6.0.0\n# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\tOrthoDB url\tDescription\n"
        (run / "full_table.tsv").write_text(hdr + "\n".join(self.ft) + "\n")
        cfg = {"prefix": "t", "organism": "Testus", "samples": [
            {"name": "S1", "isolate": "S1", "fasta": str(fa), "busco_dir": str(self.root / "busco")}]}
        cfg.update(extra_cfg or {})
        cp = self.root / "config.json"
        cp.write_text(json.dumps(cfg))
        return cp

    def seq(self, contig):
        return "".join(self.contigs[contig])


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("bp")
    ds = Dataset(root)
    exp = {}

    # 1. plus strand, single exon, complete, with stop codon: 101..193
    g1 = orf(29, 11)
    ds.place("chrA", 101, g1)
    ds.model("b1", "Complete", "chrA", "+", [(101, 193, 0)], stop=(191, 193))
    ds.row("b1", "Complete", "chrA", 100, 190, "+")
    exp["b1"] = str(Seq(g1).translate())[:-1]

    # 2. minus strand, two exons. transcript = exon1(50) + exon2(76); intron 30 nt.
    g2 = orf(40, 22)  # 126 nt
    e1, e2 = g2[:50], g2[50:]
    ds.place("chrA", 400, rc(e2) + "GT" * 15 + rc(e1))  # exon2 400..475, intron, exon1 506..555
    ds.model("b2", "Complete", "chrA", "-", [(506, 555, 0), (400, 475, 1)], stop=(400, 402))
    ds.row("b2", "Complete", "chrA", 402, 555, "-")
    exp["b2"] = str(Seq(g2).translate())[:-1]

    # 3. duplicated: two copies of the same BUSCO on the same contig (+ one on chrB)
    g3 = orf(20, 33)  # 66 nt
    ds.place("chrA", 700, g3)
    ds.place("chrA", 900, g3)
    ds.place("chrB", 300, g3)
    for s in (700, 900):
        ds.model("b3", "Duplicated", "chrA", "+", [(s, s + 65, 0)], stop=(s + 63, s + 65))
        ds.row("b3", "Duplicated", "chrA", s - 1, s + 62, "+")
    ds.model("b3", "Duplicated", "chrB", "+", [(300, 365, 0)], stop=(363, 365))
    ds.row("b3", "Duplicated", "chrB", 299, 362, "+")

    # 4. fragmented, no stop codon, + strand on chrB at 50..109 (60 nt = 20 codons)
    g4 = orf(30, 44)[:60]
    ds.place("chrB", 50, g4)
    ds.model("b4", "Fragmented", "chrB", "+", [(50, 109, 0)])
    ds.row("b4", "Fragmented", "chrB", 49, 109, "+")

    # 5. phase 1 single exon: one leading base before the first codon. 1201..1264
    g5 = "C" + orf(19, 55)  # 1 + 63 = 64 nt
    ds.place("chrA", 1201, g5)
    ds.model("b5", "Complete", "chrA", "+", [(1201, 1264, 1)], stop=(1262, 1264))
    ds.row("b5", "Complete", "chrA", 1200, 1261, "+")
    exp["b5"] = str(Seq(g5[1:]).translate())[:-1]

    # 6. frame-inconsistent two-segment model (second phase should be 2, given as 0)
    g6 = orf(30, 66)
    ds.place("chrA", 1401, g6[:40])
    ds.place("chrA", 1501, g6[40:])
    ds.model("b6", "Complete", "chrA", "+", [(1401, 1440, 0), (1501, 1501 + len(g6) - 41, 0)])
    ds.row("b6", "Complete", "chrA", 1400, 1501 + len(g6) - 44, "+")

    # 7. internal stop codon, single exon
    g7 = orf(30, 77, internal_stop_at=10)
    ds.place("chrA", 1701, g7)
    ds.model("b7", "Complete", "chrA", "+", [(1701, 1700 + len(g7), 0)], stop=(1698 + len(g7), 1700 + len(g7)))
    ds.row("b7", "Complete", "chrA", 1700, 1697 + len(g7), "+")

    # 8. miniprot Frameshift attribute (length still a multiple of 3)
    g8 = orf(30, 88)
    ds.place("chrA", 1901, g8)
    ds.model("b8", "Complete", "chrA", "+", [(1901, 1900 + len(g8), 0)], attrs="Frameshift=2")
    ds.row("b8", "Complete", "chrA", 1900, 1900 + len(g8), "+")

    # 9. rejections
    ds.model("r_oor", "Complete", "chrB", "+", [(560, 640, 0)])          # beyond chrB (600)
    ds.row("r_oor", "Complete", "chrB", 559, 640, "+")
    ds.model("r_ctg", "Complete", "chrZ", "+", [(10, 99, 0)])            # contig not in FASTA
    ds.row("r_ctg", "Complete", "chrZ", 9, 99, "+")
    ds.model("r_strand", "Complete", "chrA", "+", [(2101, 2190, 0)])     # strand mismatch
    ds.row("r_strand", "Complete", "chrA", 2100, 2190, "-")
    ds.row("r_nogff", "Complete", "chrA", 2300, 2390, "+")               # no GFF file
    ds.model("r_nomatch", "Complete", "chrA", "+", [(2501, 2590, 0)])    # GFF elsewhere
    ds.row("r_nomatch", "Complete", "chrA", 2800, 2890, "+")

    # 10. missing BUSCOs
    ds.row("m1", "Missing")
    ds.row("m2", "Missing")

    cfg = ds.write()
    out = root / "out"
    rc_ = bp.main(["--config", str(cfg), "--outdir", str(out)])
    return {"ds": ds, "out": out, "rc": rc_, "exp": exp, "cfg": cfg, "root": root}


def feats(built):
    return pd.read_csv(built["out"] / "t_S1.features.tsv", sep="\t").fillna("")


def issues(built):
    return pd.read_csv(built["out"] / "t_S1.annotation_issues.tsv", sep="\t").fillna("")


def gbk(built):
    return {r.id: r for r in SeqIO.parse(str(built["out"] / "t_S1.gbk"), "genbank")}


def cds_by_busco(built, busco_id):
    out = []
    for rec in gbk(built).values():
        for f in rec.features:
            if f.type == "CDS" and f"_{busco_id}_" in f.qualifiers["locus_tag"][0]:
                out.append((rec, f))
    return out


def test_runs_ok(built):
    assert built["rc"] == 0


def test_genbank_parses_and_sequences_identical(built):
    recs = gbk(built)
    assert set(recs) == {"chrA", "chrB"}
    for cid, rec in recs.items():
        assert str(rec.seq) == built["ds"].seq(cid)


def test_plus_strand_coordinates_exact(built):
    [(rec, f)] = cds_by_busco(built, "b1")
    assert (int(f.location.start), int(f.location.end), f.location.strand) == (100, 193, 1)
    assert f.qualifiers["translation"] == [built["exp"]["b1"]]


def test_phase_does_not_shift_boundaries(built):
    [(rec, f)] = cds_by_busco(built, "b5")
    assert (int(f.location.start), int(f.location.end)) == (1200, 1264)  # GFF 1201..1264
    assert f.qualifiers["codon_start"] == ["2"]
    assert f.qualifiers["translation"] == [built["exp"]["b5"]]


def test_minus_strand_multi_exon_join(built):
    [(rec, f)] = cds_by_busco(built, "b2")
    parts = [(int(p.start), int(p.end), p.strand) for p in f.location.parts]
    assert parts == [(505, 555, -1), (399, 475, -1)]  # transcription order
    assert f.qualifiers["codon_start"] == ["1"]
    assert f.qualifiers["translation"] == [built["exp"]["b2"]]
    genes = [g for g in rec.features if g.type == "gene" and g.qualifiers["locus_tag"] == f.qualifiers["locus_tag"]]
    assert len(genes) == 1
    assert (int(genes[0].location.start), int(genes[0].location.end), genes[0].location.strand) == (399, 555, -1)


def test_duplicated_copies_kept_distinct(built):
    hits = cds_by_busco(built, "b3")
    locs = sorted((rec.id, int(f.location.start)) for rec, f in hits)
    assert locs == [("chrA", 699), ("chrA", 899), ("chrB", 299)]
    tags = [f.qualifiers["locus_tag"][0] for _, f in hits]
    assert len(set(tags)) == 3


def test_fragmented_exported_with_status(built):
    f = feats(built)
    row = f[f.busco_id == "b4"].iloc[0]
    assert row.status == "Fragmented" and row.contig == "chrB"
    assert (row.gene_start, row.gene_end) == (50, 109)


def test_missing_never_placed(built):
    f = feats(built)
    assert not f.busco_id.isin(["m1", "m2"]).any()
    cs = pd.read_csv(built["out"] / "t_S1.contig_summary.tsv", sep="\t")
    assert not any("Missing" in c for c in cs.columns)


@pytest.mark.parametrize("bid", ["b6", "b7", "b8"])
def test_no_translation_for_invalid_frames(built, bid):
    [(rec, f)] = cds_by_busco(built, bid)
    assert "translation" not in f.qualifiers
    assert any(n.startswith("translation omitted") for n in f.qualifiers["note"])


def test_frame_inconsistent_flagged(built):
    i = issues(built)
    assert ((i.busco_id == "b6") & (i.code == "FRAME_INCONSISTENT")).any()


@pytest.mark.parametrize("bid,code", [("r_oor", "OUT_OF_RANGE"), ("r_ctg", "MISSING_CONTIG_IN_FASTA"),
                                      ("r_strand", "STRAND_MISMATCH"), ("r_nogff", "MISSING_GFF"),
                                      ("r_nomatch", "NO_OVERLAPPING_GFF_MODEL")])
def test_rejections_explicit_and_not_exported(built, bid, code):
    sk = pd.read_csv(built["out"] / "t_S1.skipped.tsv", sep="\t")
    assert ((sk.busco_id == bid) & (sk.code == code)).any()
    assert not feats(built).busco_id.eq(bid).any()


def test_no_clamping(built):
    for rec in gbk(built).values():
        for f in rec.features:
            assert 0 <= int(f.location.start) < int(f.location.end) <= len(rec)


def test_counts_orthologs_vs_hits(built):
    s = pd.read_csv(built["out"] / "t.busco_summary.tsv", sep="\t").iloc[0]
    # Complete ids: b1 b2 b5 b6 b7 b8 r_oor r_ctg r_strand r_nogff r_nomatch = 11
    assert s.busco_S_orthologs == 11 and s.busco_S_hit_rows == 11
    assert s.busco_D_orthologs == 1 and s.busco_D_hit_rows == 3
    assert s.busco_F_orthologs == 1 and s.busco_M_orthologs == 2
    assert s.busco_C_orthologs == 12
    assert s.exported_genes == 15 - 5 and s.rejected_hits == 5
    assert s.exported_genes + s.rejected_hits == s.busco_S_hit_rows + s.busco_D_hit_rows + s.busco_F_hit_rows


def test_one_gene_one_cds_per_hit(built):
    for rec in gbk(built).values():
        genes = {f.qualifiers["locus_tag"][0]: f for f in rec.features if f.type == "gene"}
        cdss = {f.qualifiers["locus_tag"][0]: f for f in rec.features if f.type == "CDS"}
        assert set(genes) == set(cdss)
        for lt, c in cdss.items():
            g = genes[lt]
            assert g.location.strand == c.location.strand
            assert int(g.location.start) == int(c.location.start) and int(g.location.end) == int(c.location.end)


def test_dry_run_writes_no_genbank(built, tmp_path):
    out = tmp_path / "dry"
    assert bp.main(["--config", str(built["cfg"]), "--outdir", str(out), "--dry-run"]) == 0
    assert not list(out.glob("*.gbk"))
    assert (out / "t_S1.features.tsv").exists() and (out / "t_S1.skipped.tsv").exists()


def test_frame_consistency_unit():
    rows = [{"start": 1, "end": 50, "phase": "0"}, {"start": 100, "end": 175, "phase": "1"}]
    assert bp.frame_consistency(rows)
    rows[1]["phase"] = "0"
    assert not bp.frame_consistency(rows)


def test_ambiguous_match_rejected():
    models = {"A": {"seqid": "c", "start": 10, "end": 100}, "B": {"seqid": "c", "start": 10, "end": 100}}
    rows = [{"Gene Start": 9, "Gene End": 100, "Sequence": "c"}]
    assigned, failed = bp.match_rows_to_models(rows, models)
    assert not assigned and failed[0].startswith("AMBIGUOUS_MATCH")


def test_full_table_header_and_missing_rows(tmp_path):
    p = tmp_path / "full_table.tsv"
    p.write_text("# BUSCO version is: 6.0.0\n# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\tOrthoDB url\tDescription\n"
                 "x1\tComplete\tc\t10\t99\t+\t1\t30\tu\td\nx2\tMissing\n")
    df = bp.read_full_table(p)
    assert list(df["Status"]) == ["Complete", "Missing"]
    assert int(df.iloc[0]["Gene Start"]) == 10
