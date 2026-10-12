# BuscoPainter

BuscoPainter **transfers BUSCO marker-gene annotations** onto the genome assembly that BUSCO was
run on. It takes the gene models that a [BUSCO](https://busco.ezlab.org/) genome-mode run already
produced (with miniprot), writes them onto the assembly as an annotated GenBank file, and reports
every model it flags or rejects. You can then open the file in
[Artemis](https://www.sanger.ac.uk/tool/artemis/) and see where each BUSCO gene sits, including
every copy of a duplicated BUSCO.

> **BuscoPainter is not a gene prediction tool, and its output is not a genome annotation.**
> It predicts nothing itself. It only re-expresses BUSCO's own output: the few thousand
> conserved marker genes of one BUSCO lineage. All other genes in the genome get no feature, and
> that is expected. The painted models are computational predictions made by BUSCO/miniprot. They
> are not curated, and their gene structures have not been verified experimentally. For a
> whole-genome annotation, use a dedicated gene prediction pipeline.

## What it does

For each sample, BuscoPainter:

1. reads the BUSCO `full_table.tsv` and the per-BUSCO GFF files that BUSCO writes in
   `busco_sequences/`;
2. matches each Complete, Duplicated or Fragmented row to exactly **one** GFF gene model
   (same BUSCO id, same contig, overlapping span);
3. writes one `gene` and one `CDS` feature per BUSCO hit, using the GFF coordinates unchanged.
   If BUSCO itself dropped terminal segments of a miniprot model (its `full_table.tsv` span and
   protein cover only some whole CDS rows), only BUSCO's segments are written, and the dropped
   ones are listed in the log;
4. adds a `/translation` only if the CDS translates cleanly in one reading frame;
5. records every flagged or rejected hit in tab-separated reports.

Nothing is silently dropped, moved or clamped. Each feature carries notes that trace it back to
its `full_table.tsv` line and its GFF `mRNA` record.

## What it does not do

- It does not run BUSCO and does not predict genes. It only re-expresses BUSCO's own output.
- It does not annotate genes other than BUSCO markers. An open reading frame with no feature is
  expected: most genes in a genome are not BUSCO markers. Every GenBank record's `source` feature
  says this in a `/note`.
- It does not correct gene models. Doubtful models (for example, multi-segment joins, frameshifts)
  are exported exactly as BUSCO's GFF gives them, and flagged.
- It does not mark partial genes with `<` or `>`, and it does not assert that a model has a true
  start or stop codon.
- It does not reproduce BUSCO's internal-stop-codon ("E") statistic.
- It does not check that your FASTA is the file BUSCO was run on. Make sure it is.

## Features

- Each Duplicated copy stays a separate, uniquely tagged feature (`<sample>_<buscoid>_<mRNA ID>`).
- Multi-segment CDS models become a single `join(...)` CDS, not several genes.
- GFF phase sets `/codon_start`. It never moves a boundary.
- Translations are attached only when they are clean, and checked by an independent validator.
- `--dry-run` runs every check and writes the reports, but no GenBank.
- The summary table gives BUSCO counts in two clearly labelled ways: ortholog counts (BUSCO's
  convention) and hit-row counts.
- A separate validator, `buscopainter-validate`, re-checks the output against the inputs without
  using the painter's code.

## Installation

You need Python 3.11 or later. Install into a fresh virtual environment:

```bash
git clone https://github.com/AroobAlhumaidy/BuscoPainter.git
cd BuscoPainter
python3 -m venv .venv
source .venv/bin/activate
pip install .
buscopainter --help
```

BuscoPainter is not on PyPI or Bioconda. Install it from this repository as shown above. It does
not need BUSCO or miniprot to be installed: it reads the files a finished BUSCO run has written.

To install the exact package versions used for validation, and run the tests:

```bash
pip install -r requirements-lock.txt -e ".[test]"
python -m pytest        # 65 tests on synthetic data; takes a few seconds
```

Dependencies: Biopython ≥ 1.88 and pandas ≥ 3.0. These are optional: PyYAML for YAML config
files (`pip install ".[yaml]"`) and requests for the OrthoDB lookup (`pip install ".[orthodb]"`).

## Quick start with the synthetic example

```bash
cd examples
buscopainter --config config.json --outdir my_output
buscopainter-validate --config config.json --outdir my_output
```

The example data are invented (random DNA, hand-made BUSCO-style files). The run takes under a
second. The outputs should match `examples/expected_output/` byte for byte. See
[examples/README.md](examples/README.md).

## Input requirements

For each sample you need:

| Input | Requirement |
|---|---|
| Assembly FASTA | **The same file BUSCO was run on.** Sequence IDs (first word of each header) must be unique and must match the `Sequence` column of `full_table.tsv`. |
| BUSCO output | A genome-mode run directory, `run_<lineage>/`, containing `full_table.tsv` and `busco_sequences/{single_copy,multi_copy,fragmented}_busco_sequences/<id>.gff`. You can give the `run_*` directory or its parent. |
| Gene predictor | **miniprot** (BUSCO's default in genome mode). Other predictors have not been tested. |

**Paths must not contain spaces.** BuscoPainter stops with an error if they do. Relative paths in
the config are resolved from the directory where you run the command.

Tested with BUSCO 6.0.0, miniprot 0.18-r281 and the `euglenozoa_odb12` lineage. Other BUSCO 5/6
versions and lineages use the same file layout, but they have not been tested. Details:
[docs/input_formats.md](docs/input_formats.md).

## Usage

Write a JSON config:

```json
{
  "prefix": "myproject",
  "organism": "Genus species",
  "samples": [
    {"name": "S1", "isolate": "S1", "fasta": "assemblies/S1.fasta", "busco_dir": "busco/S1"}
  ]
}
```

Then run:

```bash
buscopainter --config config.json --outdir results --dry-run   # checks only, no GenBank
buscopainter --config config.json --outdir results             # writes GenBank + reports
buscopainter-validate --config config.json --outdir results    # independent check
```

`python -m buscopainter` works the same way. All options and config keys are in
[docs/usage.md](docs/usage.md).

**Exit code 0 does not mean the annotation is correct.** Always read the log and the issues table,
and run the validator.

## Output files

Files are written per sample to `--outdir`:

| File | Content |
|---|---|
| `<prefix>_<sample>.gbk` | GenBank file with one record per contig, in FASTA order: `source`, `gene` and `CDS` features. |
| `<prefix>_<sample>.features.tsv` | One row per exported hit: coordinates, segments, frame, translation status and flags. |
| `<prefix>_<sample>.annotation_issues.tsv` | Every `REJECTED`, `WARNING` or `INFO` event. |
| `<prefix>_<sample>.skipped.tsv` | The `REJECTED` hits only. If it holds just a header line, nothing was rejected. |
| `<prefix>_<sample>.contig_summary.tsv` | Per-contig hit and feature counts. |
| `<prefix>.busco_summary.tsv` | Per-sample ortholog counts, hit-row counts and exported totals. |
| `<prefix>.busco_painter.log` | The run log. |

Column definitions and issue codes: [docs/output_formats.md](docs/output_formats.md).

## Viewing in Artemis

1. Install Artemis from the Wellcome Sanger Institute (it needs Java).
2. **Artemis opens only the first record of a multi-record GenBank file.** `<prefix>_<sample>.gbk`
   has one record per contig, so Artemis shows just the first contig. Run with `--per-contig` (or
   set `"export_per_contig": true`) and open the files in `<prefix>_<sample>_per_contig/`, one
   contig per file. The multi-record file is still valid GenBank (Biopython reads every record).
3. To find a BUSCO, search the `locus_tag` or `note` qualifiers for the BUSCO id (Artemis
   **Goto → Navigator**). Each CDS note gives the BUSCO status, the `full_table.tsv` line and the
   GFF evidence.
4. A multi-segment CDS (`join(...)`) is shown as linked segments. If a CDS has no `/translation`,
   its note says `translation omitted: <reason>`.

Artemis 18.2.0's own file reader was run on every per-contig file of the validation data: each
file loads the whole contig with the same features as Biopython reads.

## Validation

Release 0.1.0 (internal version `T7.1-fixed.3`) was validated on four unpublished hifiasm genome
assemblies, against BUSCO 6.0.0 (`euglenozoa_odb12`, n = 3,611, miniprot 0.18-r281). These data,
and the detailed per-sample results, are **not** distributed with this repository.

- All four samples pass all 25 checks of the independent validator,
  `buscopainter-validate --check-miniprot`. No hit was rejected, and the counts reconcile with
  BUSCO's `short_summary.json`.
- Every attached translation equals miniprot's own protein. The one known exception is a
  miniprot/BUSCO quirk: for models without a stop codon, miniprot's protein has 1 extra
  C-terminal residue.
- A second verifier that uses no BuscoPainter or Biopython code, and Artemis 18.2.0's own file
  reader, agree with every exported feature.
- The test suite passes on Python 3.11, 3.12, 3.13 and 3.14.
- The original pre-fix script exits with code 0 on the same data but fails 7 validator checks.
  This is why exit code 0 is not taken as evidence of correctness.

The validation shows that the output faithfully represents BUSCO's evidence. It does not show that
the gene models are biologically correct. Details: [docs/validation.md](docs/validation.md).

## Known limitations

These are summarised here. Full details: [docs/limitations.md](docs/limitations.md).

- **Multi-segment models are unverified.** About 10 % of miniprot models join ≥ 2 CDS segments,
  often across kilobase gaps. In organisms with few cis-introns (for example, trypanosomatids),
  these are more likely to be frameshifts, misassemblies or spurious spliced alignments. They are
  exported as given and flagged `MULTI_SEGMENT_MODEL`.
- **Translations use the standard genetic code (table 1).** Organisms with an alternative nuclear
  code are not supported.
- **Gene ends are not asserted.** Many miniprot models lack an ATG start or a stop codon.
- **Only miniprot GFFs from BUSCO genome mode** have been tested. MetaEuk, AUGUSTUS, protein mode
  and transcriptome mode have not been tested.
- **Product names** are copied verbatim from the BUSCO/OrthoDB description. They are not curated
  functional annotations.
- **Long contig IDs** (> 16 characters) give a non-standard GenBank LOCUS line.

## Citation

If you use BuscoPainter, please cite this repository (see [CITATION.cff](CITATION.cff); GitHub
shows a "Cite this repository" button). Please also cite BUSCO and miniprot, which produce the
gene models:

- Manni M., Berkeley M.R., Seppey M., Simão F.A., Zdobnov E.M. (2021). BUSCO update: novel and
  streamlined workflows along with broader and deeper phylogenetic coverage for scoring of
  eukaryotic, prokaryotic, and viral genomes. *Molecular Biology and Evolution* 38(10):4647–4654.
- Li H. (2023). Protein-to-genome alignment with miniprot. *Bioinformatics* 39(1):btad014.

## Bugs and feature requests

Please open an issue at <https://github.com/AroobAlhumaidy/BuscoPainter/issues>. Include your
BuscoPainter, Python and BUSCO versions, the command you ran, and the relevant lines of the log and
`annotation_issues.tsv`. Please do not attach unpublished genome data. See
[CONTRIBUTING.md](CONTRIBUTING.md).

## License

BuscoPainter is released under the [MIT License](LICENSE).

BUSCO and miniprot are separate tools with their own licences. They are not included in this
repository.
