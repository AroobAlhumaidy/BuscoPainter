"""Regression tests added for T7.1-fixed.2 (review of the GenBank output in Artemis, 2026-10).

Covers: strand-aware multi-segment reconstruction with every internal phase on both strands,
duplicated copies on opposite strands of one contig, GFF models absent from full_table (never
painted), partial models (rejected, never written with a segment missing), the record-level scope
note, the per-contig export used for Artemis, and the exit code on a failed count reconciliation.

Expected proteins are built from the synthetic sequences, not from the painter's own helpers.
"""
import csv

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from buscopainter import painter as bp
from test_painter import Dataset, orf, rc


def _tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _cds(out, prefix="t_S1"):
    res = {}
    for rec in SeqIO.parse(str(out / f"{prefix}.gbk"), "genbank"):
        for f in rec.features:
            if f.type == "CDS":
                res.setdefault(f.qualifiers["locus_tag"][0].split("_")[1], []).append((rec, f))
    return res


def _split(gene, cuts):
    """Split a coding sequence at the given offsets -> list of segments (transcription order)."""
    bounds = [0, *cuts, len(gene)]
    return [gene[a:b] for a, b in zip(bounds, bounds[1:])]


def _phases(segs):
    """GFF3 phase of each segment given the segments in transcription order (first phase 0)."""
    out, done = [], 0
    for s in segs:
        out.append((3 - done % 3) % 3)
        done += len(s)
    return out


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    root = tmp_path_factory.mktemp("reg")
    ds = Dataset(root)
    exp = {}

    # p3: + strand, 3 segments, internal phases 2 and 1 (cuts at 31 and 65 nt)
    g = orf(40, 101)
    segs = _split(g, [31, 65])
    ph = _phases(segs)
    assert ph == [0, 2, 1]
    starts = [1001, 1101, 1201]
    for s, seq in zip(starts, segs):
        ds.place("chrA", s, seq)
    ds.model("p3", "Complete", "chrA", "+",
             [(s, s + len(seq) - 1, p) for s, seq, p in zip(starts, segs, ph)])
    ds.row("p3", "Complete", "chrA", 1000, starts[-1] + len(segs[-1]) - 4, "+")
    exp["p3"] = str(Seq(g).translate()).rstrip("*")

    # m3: - strand, 3 segments, internal phases 1 and 2 (cuts at 32 and 64 nt)
    g = orf(40, 202)
    segs = _split(g, [32, 64])
    ph = _phases(segs)
    assert ph == [0, 1, 2]
    # transcription order runs from high to low coordinates
    ends = [1690, 1590, 1490]
    rows = []
    for e, seq, p in zip(ends, segs, ph):
        ds.place("chrA", e - len(seq) + 1, rc(seq))
        rows.append((e - len(seq) + 1, e, p))
    ds.model("m3", "Complete", "chrA", "-", rows)
    ds.row("m3", "Complete", "chrA", rows[-1][0] + 2, rows[0][1], "-")
    exp["m3"] = str(Seq(g).translate()).rstrip("*")

    # dup: Duplicated copies on opposite strands of the same contig
    g = orf(25, 303)
    ds.place("chrA", 2001, g)
    ds.place("chrA", 2201, rc(g))
    L = len(g)
    ds.model("dup", "Duplicated", "chrA", "+", [(2001, 2000 + L, 0)])
    ds.row("dup", "Duplicated", "chrA", 2000, 2000 + L - 3, "+")
    ds.model("dup", "Duplicated", "chrA", "-", [(2201, 2200 + L, 0)])
    ds.row("dup", "Duplicated", "chrA", 2203, 2200 + L, "-")
    exp["dup"] = str(Seq(g).translate()).rstrip("*")

    # extra: the GFF holds a second model that BUSCO did not report -> must not be painted
    g = orf(20, 404)
    ds.place("chrA", 2401, g)
    ds.place("chrB", 101, g)
    ds.model("extra", "Complete", "chrA", "+", [(2401, 2400 + len(g), 0)])
    ds.row("extra", "Complete", "chrA", 2400, 2400 + len(g) - 3, "+")
    ds.model("extra", "Complete", "chrB", "+", [(101, 100 + len(g), 0)])

    # partial: one CDS row has an invalid phase -> whole model rejected, not written with 1 segment
    ds.model("partial", "Complete", "chrA", "+", [(2601, 2650, 0), (2701, 2760, 1)])
    ds.row("partial", "Complete", "chrA", 2600, 2757, "+")
    ds.gff[("single_copy_busco_sequences", "partial")][-1] = \
        ds.gff[("single_copy_busco_sequences", "partial")][-1].replace("\t1\tParent=", "\t.\tParent=")

    cfg = ds.write()
    out = root / "out"
    rc_ = bp.main(["--config", str(cfg), "--outdir", str(out), "--per-contig"])
    return {"ds": ds, "out": out, "rc": rc_, "exp": exp, "cfg": cfg, "root": root}


def test_runs(run):
    assert run["rc"] == 0


@pytest.mark.parametrize("bid", ["p3", "m3"])
def test_three_segment_models_both_strands(run, bid):
    [(rec, f)] = _cds(run["out"])[bid]
    m = [l for l in run["ds"].gff[("single_copy_busco_sequences", bid)] if "\tCDS\t" in l]
    gff_rows = [(int(l.split("\t")[3]), int(l.split("\t")[4])) for l in m]
    # coordinates exactly as in the GFF, in transcription order, phases never move a boundary
    assert [(int(p.start) + 1, int(p.end)) for p in f.location.parts] == gff_rows
    assert {p.strand for p in f.location.parts} == {1 if bid == "p3" else -1}
    assert f.qualifiers["codon_start"] == ["1"]
    assert f.qualifiers["translation"] == [run["exp"][bid]]


def test_duplicates_on_opposite_strands(run):
    hits = sorted(_cds(run["out"])["dup"], key=lambda x: int(x[1].location.start))
    assert [(int(f.location.start), f.location.strand) for _, f in hits] == [(2000, 1), (2200, -1)]
    assert [f.qualifiers["translation"] for _, f in hits] == [[run["exp"]["dup"]]] * 2
    assert len({f.qualifiers["locus_tag"][0] for _, f in hits}) == 2


def test_gff_model_not_in_full_table_is_not_painted(run):
    hits = _cds(run["out"])["extra"]
    assert [(rec.id, int(f.location.start)) for rec, f in hits] == [("chrA", 2400)]
    info = [i for i in _tsv(run["out"] / "t_S1.annotation_issues.tsv")
            if i["busco_id"] == "extra" and i["code"] == "GFF_MODEL_NOT_IN_FULL_TABLE"]
    assert len(info) == 1 and info[0]["severity"] == "INFO" and "chrB:101" in info[0]["detail"]


def test_partial_model_rejected_not_written(run):
    assert "partial" not in _cds(run["out"])
    sk = _tsv(run["out"] / "t_S1.skipped.tsv")
    [r] = [r for r in sk if r["busco_id"] == "partial"]
    assert r["code"] == "INVALID_GFF_MODEL" and "phase '.'" in r["detail"]


def test_counts_reconcile_with_rejection(run):
    s = _tsv(run["out"] / "t.busco_summary.tsv")[0]
    assert s["reconciled"] == "True"
    assert int(s["exported_genes"]) + int(s["rejected_hits"]) == \
        int(s["busco_S_hit_rows"]) + int(s["busco_D_hit_rows"]) + int(s["busco_F_hit_rows"])


def test_scope_note_on_every_record(run):
    for rec in SeqIO.parse(str(run["out"] / "t_S1.gbk"), "genbank"):
        [src] = [f for f in rec.features if f.type == "source"]
        note = " ".join(src.qualifiers["note"])
        assert "only BUSCO test_odb12 marker genes are annotated" in note
        assert "not a whole-genome gene annotation" in note


def test_per_contig_files_single_record_and_identical(run):
    d = run["out"] / "t_S1_per_contig"
    files = sorted(d.glob("*.gbk"))
    assert [p.name for p in files] == ["t_S1_contig-chrA.gbk", "t_S1_contig-chrB.gbk"]
    multi = {r.id: r for r in SeqIO.parse(str(run["out"] / "t_S1.gbk"), "genbank")}
    for p in files:
        [rec] = list(SeqIO.parse(str(p), "genbank"))
        m = multi[rec.id]
        assert str(rec.seq) == str(m.seq)
        assert [(f.type, str(f.location), f.qualifiers) for f in rec.features] == \
            [(f.type, str(f.location), f.qualifiers) for f in m.features]


def test_per_contig_not_written_by_default(tmp_path, run):
    out = tmp_path / "o"
    assert bp.main(["--config", str(run["cfg"]), "--outdir", str(out)]) == 0
    assert not (out / "t_S1_per_contig").exists()
    assert "Artemis shows only the first record" in (out / "t.busco_painter.log").read_text()


def test_safe_file_part_and_collision(tmp_path):
    assert bp.safe_file_part("ctg/1|x y") == "ctg_1_x_y"
    assert bp.safe_file_part("../etc") == "_etc"
    from Bio.SeqRecord import SeqRecord
    recs = {k: SeqRecord(Seq("ACGT"), id=k, annotations={"molecule_type": "DNA"}) for k in ("a/b", "a|b")}
    with pytest.raises(bp.PainterError, match="same file name"):
        bp.write_per_contig_genbanks(tmp_path, "p", "s", recs)


def test_failed_reconciliation_returns_3(tmp_path, run, monkeypatch):
    real = bp.add_busco_features

    def drop_one(*a, **k):
        feats, issues = real(*a, **k)
        return feats[1:], issues  # simulate a silently lost hit

    monkeypatch.setattr(bp, "add_busco_features", drop_one)
    assert bp.main(["--config", str(run["cfg"]), "--outdir", str(tmp_path / "o")]) == 3
