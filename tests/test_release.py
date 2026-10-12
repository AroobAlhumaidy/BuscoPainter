"""Release tests: packaging, CLI, the synthetic example, and extra edge cases.

Expected values come from the synthetic sequences (examples/make_example_data.py and the Dataset
helper in test_painter.py), not from the painter's own helpers.
"""
import csv
import hashlib
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from buscopainter import painter, validate
from test_painter import Dataset, orf

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "examples"
EXPECTED = EXAMPLES / "expected_output"

# SHA-256 of the validated reference files: painter.py = T7.1-fixed.3 (see CHANGELOG 0.1.0 and
# docs/validation.md); validate.py = validate_busco_painter.py + the T7.1-fixed.3
# BUSCO-trimmed-model checks.
# Changing painter.py changes the scientific reference: update this hash only together with a
# CHANGELOG entry and a fresh validation on real data (see docs/validation.md).
REFERENCE_SHA256 = {
    "painter.py": "494be47bd68583f9bf43e1a5c390d3734966035376b1ed344d42f982bd42d235",
    "validate.py": "d357e99fe74e13d54ac59220db54b56d23428e4c3d6110f181201a02418a151f",
}


def _tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _load_generator():
    spec = importlib.util.spec_from_file_location("make_example_data", EXAMPLES / "make_example_data.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------- packaging ---
@pytest.mark.parametrize("module", [painter, validate])
def test_reference_implementation_unchanged(module):
    name = Path(module.__file__).name
    assert hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() == REFERENCE_SHA256[name]


@pytest.mark.parametrize("cmd", [[sys.executable, "-m", "buscopainter", "--help"],
                                 [sys.executable, "-m", "buscopainter.validate", "--help"]])
def test_cli_help(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0
    assert "--config" in r.stdout


def test_missing_config_returns_error(tmp_path):
    assert painter.main(["--config", str(tmp_path / "nope.json"), "--outdir", str(tmp_path / "o")]) == 2


# ------------------------------------------------------- synthetic example ---
@pytest.fixture(scope="module")
def example_run(tmp_path_factory):
    """Copy examples/ to a temp dir and run from there (config paths are relative)."""
    work = tmp_path_factory.mktemp("example")
    shutil.copy(EXAMPLES / "config.json", work)
    shutil.copytree(EXAMPLES / "data", work / "data")
    mp = pytest.MonkeyPatch()
    mp.chdir(work)
    try:
        rc = painter.main(["--config", "config.json", "--outdir", "out"])
        rc_dry = painter.main(["--config", "config.json", "--outdir", "dry", "--dry-run"])
        rc_val = validate.main(["--config", "config.json", "--outdir", "out"])
    finally:
        mp.undo()
    return {"work": work, "out": work / "out", "dry": work / "dry", "rc": rc, "rc_dry": rc_dry, "rc_val": rc_val}


def test_example_data_is_reproducible(tmp_path, monkeypatch):
    gen = _load_generator()
    monkeypatch.setattr(gen, "DATA", tmp_path / "data")
    gen.main()
    for f in sorted((EXAMPLES / "data").rglob("*")):
        if f.is_file():
            assert (tmp_path / "data" / f.relative_to(EXAMPLES / "data")).read_bytes() == f.read_bytes(), f


def test_example_runs_and_validates(example_run):
    assert example_run["rc"] == 0
    assert example_run["rc_val"] == 0


@pytest.mark.parametrize("name", sorted(p.name for p in EXPECTED.iterdir()))
def test_example_matches_expected_output(example_run, name):
    assert (example_run["out"] / name).read_bytes() == (EXPECTED / name).read_bytes()


def _example_cds(out):
    cds = {}
    for rec in SeqIO.parse(str(out / "example_demo.gbk"), "genbank"):
        for f in rec.features:
            if f.type == "CDS":
                cds[f.qualifiers["locus_tag"][0]] = (rec, f)
    return cds


def test_example_translations_equal_designed_proteins(example_run):
    gen = _load_generator()
    designed = {"ex1at0": gen.orf(60, 1), "ex2at0": gen.orf(80, 2), "ex3at0": gen.orf(40, 3),
                "ex4at0": gen.orf(45, 4), "ex5at0": gen.orf(70, 5)[:150]}
    want = {k: str(Seq(v).translate()).rstrip("*") for k, v in designed.items()}
    seen = set()
    for lt, (rec, f) in _example_cds(example_run["out"]).items():
        bid = lt.split("_")[1]
        if bid in want:
            assert f.qualifiers["translation"] == [want[bid]], lt
            seen.add(bid)
        else:  # ex6 (frameshift), ex7 (frame-inconsistent)
            assert "translation" not in f.qualifiers, lt
    assert seen == set(want)


def test_example_reverse_strand_and_phase(example_run):
    cds = {lt.split("_")[1]: f for lt, (_, f) in _example_cds(example_run["out"]).items()}
    # ex2: GFF CDS 1037..1086 (phase 0) then 801..996 (phase 1), minus strand, transcription order
    assert [(int(p.start), int(p.end), p.strand) for p in cds["ex2at0"].location.parts] == \
        [(1036, 1086, -1), (800, 996, -1)]
    # ex3: GFF 1501..1627 phase 1 -> boundaries unchanged, codon_start 2
    assert (int(cds["ex3at0"].location.start), int(cds["ex3at0"].location.end)) == (1500, 1627)
    assert cds["ex3at0"].qualifiers["codon_start"] == ["2"]


def test_example_duplicates_same_contig_distinct(example_run):
    copies = sorted((rec.id, int(f.location.start), f.location.strand)
                    for lt, (rec, f) in _example_cds(example_run["out"]).items() if "_ex4at0_" in lt)
    assert copies == [("contig_1", 2200, 1), ("contig_1", 2800, 1), ("contig_2", 300, -1)]


def test_example_missing_rows_counted_not_placed(example_run):
    tags = set(_example_cds(example_run["out"]))
    assert not any("_ex8at0_" in t or "_ex9at0_" in t for t in tags)
    s = _tsv(example_run["out"] / "example.busco_summary.tsv")[0]
    assert s["busco_M_orthologs"] == "2"


def test_example_locus_tags_unique(example_run):
    for ftype in ("gene", "CDS"):
        tags = [f.qualifiers["locus_tag"][0]
                for rec in SeqIO.parse(str(example_run["out"] / "example_demo.gbk"), "genbank")
                for f in rec.features if f.type == ftype]
        assert len(tags) == len(set(tags)) == 9


def test_example_output_files_consistent(example_run):
    out = example_run["out"]
    feats = _tsv(out / "example_demo.features.tsv")
    cds = _example_cds(out)
    assert {r["locus_tag"] for r in feats} == set(cds)
    for r in feats:
        rec, f = cds[r["locus_tag"]]
        assert r["contig"] == rec.id
        assert r["cds_segments"] == str(len(f.location.parts))
        assert (r["translation_attached"] == "True") == ("translation" in f.qualifiers)
    summ = _tsv(out / "example.busco_summary.tsv")[0]
    total = [r for r in _tsv(out / "example_demo.contig_summary.tsv") if r["contig"] == "TOTAL"][0]
    assert summ["exported_genes"] == summ["exported_CDS"] == total["exported_genes"] == str(len(feats))
    assert summ["cds_segments"] == total["cds_segments"] == str(sum(int(r["cds_segments"]) for r in feats))
    issues = _tsv(out / "example_demo.annotation_issues.tsv")
    skipped = _tsv(out / "example_demo.skipped.tsv")
    assert skipped == [i for i in issues if i["severity"] == "REJECTED"]


def test_dry_run_same_reports_no_genbank(example_run):
    assert example_run["rc_dry"] == 0
    assert not list(example_run["dry"].glob("*.gbk"))
    for name in ("example_demo.features.tsv", "example_demo.annotation_issues.tsv",
                 "example_demo.skipped.tsv", "example_demo.contig_summary.tsv"):
        assert (example_run["dry"] / name).read_bytes() == (example_run["out"] / name).read_bytes()
    full, dry = (_tsv(example_run[k] / "example.busco_summary.tsv")[0] for k in ("out", "dry"))
    assert dry.pop("dry_run") == "True" and full.pop("dry_run") == "False"
    assert dry == full


# ----------------------------------------------------- contig boundaries ---
def test_contig_boundary_exact_end_accepted_one_past_rejected(tmp_path):
    ds = Dataset(tmp_path)                      # chrB is 600 bp
    g = orf(28, 5)                              # 90 nt
    ds.place("chrB", 511, g)
    ds.model("edge_end", "Complete", "chrB", "+", [(511, 600, 0)], stop=(598, 600))
    ds.row("edge_end", "Complete", "chrB", 510, 597, "+")
    ds.place("chrA", 1, str(Seq(g).reverse_complement()))
    ds.model("edge_start", "Complete", "chrA", "-", [(1, 90, 0)], stop=(1, 3))
    ds.row("edge_start", "Complete", "chrA", 3, 90, "-")
    ds.model("past_end", "Complete", "chrB", "+", [(512, 601, 0)])
    ds.row("past_end", "Complete", "chrB", 511, 601, "+")
    cfg = ds.write()
    out = tmp_path / "out"
    assert painter.main(["--config", str(cfg), "--outdir", str(out)]) == 0
    feats = {r["busco_id"]: r for r in _tsv(out / "t_S1.features.tsv")}
    assert (feats["edge_end"]["gene_start"], feats["edge_end"]["gene_end"]) == ("511", "600")
    assert (feats["edge_start"]["gene_start"], feats["edge_start"]["gene_end"]) == ("1", "90")
    want = str(Seq(g).translate()).rstrip("*")
    bid_of = {r["locus_tag"]: bid for bid, r in feats.items()}
    cds = {bid_of[f.qualifiers["locus_tag"][0]]: f
           for rec in SeqIO.parse(str(out / "t_S1.gbk"), "genbank") for f in rec.features if f.type == "CDS"}
    assert cds["edge_end"].qualifiers["translation"] == [want]
    assert cds["edge_start"].qualifiers["translation"] == [want]
    sk = _tsv(out / "t_S1.skipped.tsv")
    assert [(r["busco_id"], r["code"]) for r in sk] == [("past_end", "OUT_OF_RANGE")]
    assert validate.main(["--config", str(cfg), "--outdir", str(out)]) == 0
