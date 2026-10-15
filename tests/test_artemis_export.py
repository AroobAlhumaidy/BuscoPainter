"""Tests for the whole-assembly Artemis export (indexed FASTA + indexed GFF3).

Expected values come from the designed synthetic example (examples/make_example_data.py and
examples/expected_output/), worked out by hand, not from the exporter's own helpers.
"""
import gzip
import shutil
from pathlib import Path

import pytest

from buscopainter import artemis_export, painter, validate_artemis

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "examples"
HAS_TABIX = shutil.which("tabix") is not None
needs_tabix = pytest.mark.skipif(not HAS_TABIX, reason="tabix (htslib) not installed")

# locus_tag -> (contig, strand, gene start, gene end, [(CDS start, end, phase) in transcription order])
DESIGNED = {
    "demo_ex1at0_MP000001": ("contig_1", "+", 201, 386, [(201, 386, 0)]),
    # minus-strand join: 1037..1086 (50 bp) is read first; 50 % 3 = 2 -> next segment phase 1
    "demo_ex2at0_MP000002": ("contig_1", "-", 801, 1086, [(1037, 1086, 0), (801, 996, 1)]),
    "demo_ex3at0_MP000003": ("contig_1", "+", 1501, 1627, [(1501, 1627, 1)]),       # codon_start 2
    "demo_ex4at0_MP000004": ("contig_1", "+", 2201, 2341, [(2201, 2341, 0)]),
    "demo_ex4at0_MP000005": ("contig_1", "+", 2801, 2941, [(2801, 2941, 0)]),
    "demo_ex4at0_MP000006": ("contig_2", "-", 301, 441, [(301, 441, 0)]),
    "demo_ex5at0_MP000007": ("contig_2", "+", 1001, 1150, [(1001, 1150, 0)]),
    "demo_ex6at0_MP000008": ("contig_1", "+", 3601, 3756, [(3601, 3756, 0)]),
    # 4201..4240 is 40 bp; 40 % 3 = 1 -> next segment phase 2
    "demo_ex7at0_MP000009": ("contig_1", "+", 4201, 4546, [(4201, 4240, 0), (4401, 4546, 2)]),
}


@pytest.fixture()
def run(tmp_path, monkeypatch):
    shutil.copy(EXAMPLES / "config.json", tmp_path)
    shutil.copytree(EXAMPLES / "data", tmp_path / "data")
    monkeypatch.chdir(tmp_path)
    assert painter.main(["--config", "config.json", "--outdir", "out"]) == 0
    rc = artemis_export.main(["--config", "config.json", "--outdir", "out"])
    return {"dir": tmp_path, "exp": tmp_path / "out" / "example_demo_artemis", "rc": rc}


def _rows(gff):
    with gzip.open(gff, "rt") as fh:
        return [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]


def _validate(run):
    return validate_artemis.main(["--config", "config.json", "--outdir", "out"])


def test_export_files_and_fasta_unchanged(run):
    assert run["rc"] == 0
    exp = run["exp"]
    src = run["dir"] / "data" / "genome.fasta"
    assert (exp / "genome.fasta").read_bytes() == src.read_bytes()
    assert (exp / "genome.fasta.fai").exists()
    assert (exp / "example_demo.busco.gff3.gz").exists()
    assert not (run["dir"] / "data" / "genome.fasta.fai").exists()   # input folder untouched


def test_gff3_models_equal_designed(run):
    rows = _rows(run["exp"] / "example_demo.busco.gff3.gz")
    got = {}
    for f in rows:
        lt = dict(kv.split("=", 1) for kv in f[8].split(";"))["locus_tag"]
        m = got.setdefault(lt, {"CDS": []})
        if f[2] == "CDS":
            m["CDS"].append((int(f[3]), int(f[4]), int(f[7])))
        else:
            m[f[2]] = (f[0], f[6], int(f[3]), int(f[4]))
    assert set(got) == set(DESIGNED)
    for lt, (contig, strand, gs, ge, cds) in DESIGNED.items():
        assert got[lt]["gene"] == got[lt]["mRNA"] == (contig, strand, gs, ge), lt
        assert sorted(got[lt]["CDS"]) == sorted(cds), lt


def test_sequence_regions_and_sorting(run):
    with gzip.open(run["exp"] / "example_demo.busco.gff3.gz", "rt") as fh:
        text = fh.read()
    assert text.startswith("##gff-version 3\n")
    fai = (run["exp"] / "genome.fasta.fai").read_text().splitlines()
    for line in fai:
        name, length = line.split("\t")[:2]
        assert f"##sequence-region {name} 1 {length}\n" in text
    rows = _rows(run["exp"] / "example_demo.busco.gff3.gz")
    order = [l.split("\t")[0] for l in fai]
    keys = [(order.index(f[0]), int(f[3])) for f in rows]
    assert keys == sorted(keys)


def test_notes_are_escaped_and_round_trip(run):
    rows = _rows(run["exp"] / "example_demo.busco.gff3.gz")
    for f in rows:
        for kv in f[8].split(";"):
            k, v = kv.split("=", 1)
            assert "=" not in v and ";" not in v
    note = [f[8] for f in rows if f[2] == "CDS" and "demo_ex1at0" in f[8] and "Note=" in f[8]][0]
    from urllib.parse import unquote
    first = unquote(note.split("Note=")[1].split(",")[0])
    assert first.startswith("BUSCO ex1at0 (Complete), lineage example_odb12; full_table line 5:")


@needs_tabix
def test_validator_passes(run):
    assert (run["exp"] / "example_demo.busco.gff3.gz.tbi").exists()
    assert _validate(run) == 0


@needs_tabix
@pytest.mark.parametrize("mutation", ["cds_coord", "strand", "drop_model", "fasta_base", "no_tbi", "phase"])
def test_validator_catches_mutations(run, mutation):
    exp = run["exp"]
    gff = exp / "example_demo.busco.gff3.gz"
    with gzip.open(gff, "rt") as fh:
        lines = fh.read().splitlines(keepends=True)
    if mutation == "fasta_base":
        fa = exp / "genome.fasta"
        b = bytearray(fa.read_bytes())
        i = b.index(b"\n") + 10
        b[i] = ord("A") if b[i] != ord("A") else ord("C")
        fa.write_bytes(bytes(b))
    elif mutation == "no_tbi":
        Path(str(gff) + ".tbi").unlink()
    else:
        for i, l in enumerate(lines):
            f = l.split("\t")
            if mutation == "cds_coord" and len(f) == 9 and f[2] == "CDS" and "ex7at0" in f[8] and f[3] == "4401":
                f[3] = "4404"
            elif mutation == "strand" and len(f) == 9 and "ex1at0" in f[8]:
                f[6] = "-"
            elif mutation == "phase" and len(f) == 9 and f[2] == "CDS" and f[3] == "801":
                f[7] = "0"
            elif mutation == "drop_model" and len(f) == 9 and "ex5at0" in f[8]:
                lines[i] = ""
                continue
            lines[i] = "\t".join(f)
        from Bio import bgzf
        with bgzf.BgzfWriter(str(gff), "wb") as w:
            w.write("".join(lines).encode())
        artemis_export.subprocess.run([shutil.which("tabix"), "-f", "-p", "gff", str(gff)], check=True)
    assert _validate(run) == 1


def test_fasta_index_matches_samtools_rules(tmp_path):
    fa = tmp_path / "x.fasta"
    fa.write_bytes(b">a desc\nACGTA\nCGTAC\nGG\n>b\nTTTT\n")
    assert artemis_export.fasta_index(fa) == [("a", 12, 8, 5, 6), ("b", 4, 26, 4, 5)]
    bad = tmp_path / "bad.fasta"
    bad.write_bytes(b">a\nACG\nACGTA\n")
    with pytest.raises(artemis_export.ExportError):
        artemis_export.fasta_index(bad)
    bad.write_bytes(b">a\nACGTA\nAC\nACGTA\n")
    with pytest.raises(artemis_export.ExportError):
        artemis_export.fasta_index(bad)


def test_export_refuses_a_different_fasta(run, tmp_path):
    other = tmp_path / "other.fasta"
    text = (run["dir"] / "data" / "genome.fasta").read_text().splitlines()
    text[1] = ("T" if text[1][0] != "T" else "G") + text[1][1:]
    other.write_text("\n".join(text) + "\n")
    rc = artemis_export.main(["--gbk", str(run["dir"] / "out" / "example_demo.gbk"),
                              "--fasta", str(other), "--outdir", str(tmp_path / "e")])
    assert rc == 2


def test_without_tabix_errors_clearly(tmp_path, monkeypatch, capsys):
    shutil.copy(EXAMPLES / "config.json", tmp_path)
    shutil.copytree(EXAMPLES / "data", tmp_path / "data")
    monkeypatch.chdir(tmp_path)
    assert painter.main(["--config", "config.json", "--outdir", "out"]) == 0
    monkeypatch.setattr(artemis_export.shutil, "which", lambda name: None)
    assert artemis_export.main(["--config", "config.json", "--outdir", "out"]) == 2
    err = capsys.readouterr().err
    assert "tabix not found" in err
    assert "Install htslib/tabix" in err
    assert not (tmp_path / "out" / "example_demo_artemis").exists()
