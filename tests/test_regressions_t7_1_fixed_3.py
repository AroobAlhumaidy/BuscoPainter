"""Regression tests added for T7.1-fixed.3 (final verification, 2026-10-12).

BUSCO 6 (miniprot mode) sometimes drops terminal miniprot CDS segments: full_table.tsv and the
BUSCO .faa protein then cover only the remaining segments, while the per-BUSCO GFF still holds the
whole model. T7.1-fixed.1/.2 painted the whole GFF model (up to 10.8 kb of segments BUSCO had
discarded). These tests pin the corrected behaviour:
paint BUSCO's prediction, log what was dropped, and fall back to the old behaviour when the
full_table span does not fall on whole segments.
"""
import csv

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from buscopainter import painter as bp
from buscopainter import validate as bv
from test_painter import Dataset, orf, rc


def _tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _cds(out, prefix="t_S1"):
    res = {}
    for rec in SeqIO.parse(str(out / f"{prefix}.gbk"), "genbank"):
        for f in rec.features:
            if f.type in ("gene", "CDS"):
                res.setdefault((f.qualifiers["locus_tag"][0].split("_")[1], f.type), []).append(f)
    return res


def _tag_cds_row(ds, bid, start, attr):
    key = ("single_copy_busco_sequences", bid)
    ds.gff[key] = [l.replace("\tParent=", f"\t{attr};Parent=") if f"\tCDS\t{start}\t" in l else l
                   for l in ds.gff[key]]


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    root = tmp_path_factory.mktemp("reg3")
    ds = Dataset(root)
    exp = {}

    # tp: + strand; spurious 5' segment far upstream, BUSCO span = the real segment only
    g = orf(40, 11)
    ds.place("chrA", 501, g)
    ds.model("tp", "Complete", "chrA", "+", [(101, 130, 0), (501, 500 + len(g), 0)])
    ds.row("tp", "Complete", "chrA", 501, 500 + len(g), "+")
    exp["tp"] = (501, 500 + len(g), 1, str(Seq(g).translate()).rstrip("*"))

    # tpp: + strand; the kept segment starts with phase 2 (2 bases finish a codon begun upstream)
    body = orf(30, 12)
    seg = "GC" + body
    ds.place("chrA", 701, seg)
    ds.model("tpp", "Complete", "chrA", "+", [(201, 231, 0), (701, 700 + len(seg), 2)])
    ds.row("tpp", "Complete", "chrA", 701, 700 + len(seg), "+")
    exp["tpp"] = (701, 700 + len(seg), 3, str(Seq(body).translate()).rstrip("*"))

    # tm: - strand; spurious 3' segment (lower coordinates) carries miniprot's Frameshift tag
    g = orf(40, 13)
    ds.place("chrA", 1201, rc(g))
    ds.model("tm", "Complete", "chrA", "-", [(1201, 1200 + len(g), 0), (1001, 1030, 0)],
             attrs="Frameshift=1")
    _tag_cds_row(ds, "tm", 1001, "Frameshift=1")
    ds.row("tm", "Complete", "chrA", 1201, 1200 + len(g), "-")
    exp["tm"] = (1201, 1200 + len(g), 1, str(Seq(g).translate()).rstrip("*"))

    # tcut: full_table span cuts through a segment -> not trimmed; whole model, no translation
    g = orf(40, 14)
    ds.place("chrA", 1601, g)
    ds.model("tcut", "Complete", "chrA", "+", [(1501, 1530, 0), (1601, 1600 + len(g), 0)])
    ds.row("tcut", "Complete", "chrA", 1520, 1600 + len(g), "+")

    cfg = ds.write()
    out = root / "out"
    rc_ = bp.main(["--config", str(cfg), "--outdir", str(out), "--per-contig"])
    return {"ds": ds, "out": out, "rc": rc_, "exp": exp, "cfg": cfg}


def test_runs(run):
    assert run["rc"] == 0


@pytest.mark.parametrize("bid", ["tp", "tpp", "tm"])
def test_trimmed_to_busco_span(run, bid):
    s, e, cs, prot = run["exp"][bid]
    feats = _cds(run["out"])
    [cds], [gene] = feats[(bid, "CDS")], feats[(bid, "gene")]
    assert [(int(p.start) + 1, int(p.end)) for p in cds.location.parts] == [(s, e)]
    assert (int(gene.location.start) + 1, int(gene.location.end)) == (s, e)
    assert cds.qualifiers["codon_start"] == [str(cs)]
    assert cds.qualifiers["translation"] == [prot]
    assert any("were not exported" in n for n in cds.qualifiers["note"])


def test_trim_is_logged(run):
    iss = [i for i in _tsv(run["out"] / "t_S1.annotation_issues.tsv") if i["code"] == "FULLTABLE_SPAN_TRIMMED"]
    assert sorted(i["busco_id"] for i in iss) == ["tm", "tp", "tpp"]
    assert all(i["severity"] == "WARNING" for i in iss)
    assert "not exported: 101-130" in next(i for i in iss if i["busco_id"] == "tp")["detail"]
    ft = {r["busco_id"]: r for r in _tsv(run["out"] / "t_S1.features.tsv")}
    assert "FULLTABLE_SPAN_TRIMMED" in ft["tp"]["flags"]
    assert (ft["tp"]["gene_start"], ft["tp"]["cds_segments"]) == ("501", "1")


def test_span_not_on_segment_boundary_keeps_whole_model(run):
    feats = _cds(run["out"])
    [cds] = feats[("tcut", "CDS")]
    assert len(cds.location.parts) == 2
    assert "translation" not in cds.qualifiers
    codes = {i["code"] for i in _tsv(run["out"] / "t_S1.annotation_issues.tsv") if i["busco_id"] == "tcut"}
    assert "FULLTABLE_SPAN_DIFFERS_FROM_GFF" in codes and "FULLTABLE_SPAN_TRIMMED" not in codes


def test_validator_accepts_trimmed_output(run):
    assert bv.main(["--config", str(run["cfg"]), "--outdir", str(run["out"])]) == 0


def test_validator_rejects_untrimmed_model(run, tmp_path):
    """If the painter wrote the whole GFF model again, the validator must fail."""
    import shutil
    out = tmp_path / "out"
    shutil.copytree(run["out"], out)
    gbk = out / "t_S1.gbk"
    s, e, _, _ = run["exp"]["tp"]
    txt = gbk.read_text()
    assert txt.count(f"{s}..{e}") == 2  # gene + CDS
    gbk.write_text(txt.replace(f"CDS             {s}..{e}", f"CDS             join(101..130,{s}..{e})", 1))
    assert bv.main(["--config", str(run["cfg"]), "--outdir", str(out)]) == 1
