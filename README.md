# BuscoPainter

BuscoPainter transfers BUSCO marker-gene annotations onto a genome assembly.

It writes GenBank for normal use.

It can also write indexed FASTA + GFF3 for Artemis.

Version 0.1.0 - Python 3.11 or later - [MIT License](LICENSE)

## What Is BuscoPainter?

[BUSCO](https://busco.ezlab.org/) checks whether expected marker genes are present in a genome.

In genome mode, BUSCO uses miniprot to place those marker genes on the assembly.

BUSCO saves those placements in many small files.

BuscoPainter collects them.

It paints them back onto the original FASTA assembly.

## What Does It Do?

BuscoPainter writes:

- one combined GenBank file;
- optional one-contig GenBank files;
- optional Artemis files for the whole assembly;
- TSV reports for exported, flagged and rejected hits.

Each exported BUSCO hit becomes one `gene` feature and one `CDS` feature.

Duplicated BUSCOs are kept as separate copies.

## What It Does Not Do

BuscoPainter does not run BUSCO.

BuscoPainter does not predict new genes.

BuscoPainter is not a whole-genome annotation pipeline.

Only BUSCO marker genes are annotated.

Most real genes are not BUSCO markers.

So unannotated genes and ORFs are expected.

BUSCO completeness is not annotation completeness.

Always validate the output before using it.

## Quick Start

These commands install BuscoPainter and run the safe synthetic example.

The example data are invented. They contain no real genome data.

### 1. Install BuscoPainter

```bash
git clone https://github.com/AroobAlhumaidy/BuscoPainter.git
cd BuscoPainter
python3 -m venv .venv
source .venv/bin/activate
pip install .
```

### 2. Run the Example

```bash
cd examples
buscopainter --config config.json --outdir my_output
```

This writes the combined GenBank file:

```text
my_output/example_demo.gbk
```

### 3. Validate the GenBank Output

```bash
buscopainter-validate --config config.json --outdir my_output
```

The last line should be:

```text
OVERALL: PASS
```

### 4. Install `tabix` for Artemis

`tabix` creates the GFF3 index that Artemis needs.

With conda:

```bash
conda install -c bioconda htslib
```

With Homebrew on macOS:

```bash
brew install htslib
```

Check it:

```bash
tabix --version
```

### 5. Make Artemis Files

```bash
buscopainter-artemis --config config.json --outdir my_output
```

This writes:

```text
my_output/example_demo_artemis/genome.fasta
my_output/example_demo_artemis/genome.fasta.fai
my_output/example_demo_artemis/example_demo.busco.gff3.gz
my_output/example_demo_artemis/example_demo.busco.gff3.gz.tbi
```

### 6. Validate the Artemis Files

```bash
buscopainter-validate-artemis --config config.json --outdir my_output
```

The last line should be:

```text
OVERALL: PASS
```

## Run BuscoPainter On Your Own Genome

You need two inputs for each sample.

First, use the same FASTA file that BUSCO used.

Second, use the BUSCO genome-mode output folder.

The folder should contain a `run_<lineage>` directory.

That directory should contain `full_table.tsv` and `busco_sequences/`.

### 1. Make `config.json`

Edit the paths for your data.

```json
{
  "prefix": "myproject",
  "organism": "Genus species",
  "samples": [
    {
      "name": "S1",
      "isolate": "S1",
      "fasta": "assemblies/S1.fasta",
      "busco_dir": "busco/S1"
    }
  ]
}
```

Paths are read from the folder where you run the command.

Paths must not contain spaces.

### 2. Do A Dry Run

```bash
buscopainter --config config.json --outdir results --dry-run
```

Open this file:

```text
results/myproject_S1.annotation_issues.tsv
```

Check for rejected hits.

### 3. Write GenBank

```bash
buscopainter --config config.json --outdir results
```

This writes:

```text
results/myproject_S1.gbk
```

### 4. Optional: Write One GenBank Per Contig

Use this if you prefer separate Artemis files.

```bash
buscopainter --config config.json --outdir results --per-contig
```

The files go here:

```text
results/myproject_S1_per_contig/
```

### 5. Validate

```bash
buscopainter-validate --config config.json --outdir results --report results/validation.json
```

For a deeper translation check, add `--check-miniprot`.

That needs the BUSCO parent folder with `logs/`.

## Open Your Complete Genome In Artemis

[Artemis](https://www.sanger.ac.uk/tool/artemis/) is a genome browser.

Artemis needs Java.

For a whole assembly in one Artemis window, BuscoPainter uses:

- indexed multi-FASTA;
- BGZF-compressed GFF3;
- tabix `.tbi` index.

Create those files:

```bash
buscopainter-artemis --config config.json --outdir results
```

Validate them:

```bash
buscopainter-validate-artemis --config config.json --outdir results
```

Open them in Artemis:

1. Choose **File -> Open ...**.
2. Open `results/myproject_S1_artemis/S1.fasta`.
3. Click **Use index**.
4. Choose **File -> Read An Entry ...**.
5. Open `results/myproject_S1_artemis/myproject_S1.busco.gff3.gz`.
6. Click **Use index**.
7. Choose a contig from the drop-down list.

Do not use **Concatenate Sequences**.

Do not open the combined `.gbk` if you need all contigs in Artemis.

Artemis shows only the first contig of a multi-record GenBank file.

The combined `.gbk` is still useful for other tools.

More detail is in [docs/artemis.md](docs/artemis.md).

## Understand Your Results

Start with these files:

- `myproject_S1.gbk`: BUSCO marker genes on the assembly.
- `myproject_S1.features.tsv`: one row per exported BUSCO hit.
- `myproject_S1.annotation_issues.tsv`: warnings and rejected hits.
- `myproject_S1.skipped.tsv`: rejected hits only.
- `myproject.busco_summary.tsv`: counts for each sample.

Issue levels:

- `REJECTED` means the hit was not exported.
- `WARNING` means the hit was exported but needs review.
- `INFO` means a note only.

`OVERALL: PASS` means the files match BUSCO's evidence.

It does not mean the gene model is biologically correct.

Output columns and issue codes are explained in [docs/output_formats.md](docs/output_formats.md).

## Important Limitations

BuscoPainter transfers existing BUSCO marker annotations.

It does not create a complete genome annotation.

Multi-segment models are exported as BUSCO/miniprot gave them.

They may represent real introns, frameshifts, assembly errors or false alignments.

Gene starts and ends are not curated.

The standard genetic code is used for translation.

Product names come from BUSCO descriptions.

They are not curated functional annotations.

Always inspect the reports.

Always validate important results.

The full list is in [docs/limitations.md](docs/limitations.md).

## Documentation

- [Usage guide](docs/usage.md): commands, config keys and exit codes.
- [Input formats](docs/input_formats.md): required BUSCO files.
- [Output formats](docs/output_formats.md): output files and issue codes.
- [Artemis guide](docs/artemis.md): whole-assembly viewing.
- [Limitations](docs/limitations.md): scientific and technical limits.
- [Validation report](docs/validation.md): how this release was checked.
- [Synthetic example](examples/README.md): what the example tests.
- [Changelog](CHANGELOG.md): release history.

## Tested Setup

BuscoPainter has been tested with Python 3.11 to 3.14.

The validation data used BUSCO 6.0.0 in genome mode with miniprot 0.18-r281.

The whole-assembly Artemis workflow was tested with Artemis 18.2.0.

The validation assemblies are unpublished and are not included.

## Run The Tests

```bash
pip install -r requirements-lock.txt -e ".[test]"
python -m pytest
```

Some Artemis tests need `tabix`.

Without `tabix`, those tests are skipped.

## Citation And License

If you use BuscoPainter, please cite this repository.

Citation details are in [CITATION.cff](CITATION.cff).

Please also cite BUSCO and miniprot.

BuscoPainter is released under the [MIT License](LICENSE).

BUSCO and miniprot are separate tools with their own licenses.
