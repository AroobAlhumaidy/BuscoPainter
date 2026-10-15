# Validation

## Scope

Release 0.1.0 packages internal version `T7.1-fixed.3`: the reference script
`buscoPainter.T7.1.fixed.py` (`T7.1-fixed.1`) plus reporting and Artemis changes (`T7.1-fixed.2`)
and one scientific correction: models that BUSCO trimmed to some whole CDS rows are exported as
BUSCO's prediction (`T7.1-fixed.3`; see below and the CHANGELOG). The test
`test_reference_implementation_unchanged` checks the SHA-256 of `painter.py`.

The validation shows that the output faithfully represents the BUSCO/miniprot evidence. Every
feature traces to its input rows, the coordinates equal the GFF, the sequences are unchanged, the
translations are clean, and the counts reconcile. **It does not show that the gene models are
biologically correct.** No experimental or curation-based validation of gene structures has been
done.

## Real-data validation

Data: four unpublished hifiasm primary genome assemblies. The FASTA files are the files BUSCO was
run on (SHA-256 checksums identical). BUSCO 6.0.0, `euglenozoa_odb12` (n = 3,611),
miniprot 0.18-r281. The data and the per-sample results are not distributed with this repository.

Summary of the results for `T7.1-fixed.3`:

| Check | Result |
|---|---|
| `buscopainter-validate --check-miniprot` | 4/4 samples PASS, 100/100 checks |
| `buscopainter-validate --check-faa --check-miniprot` (one sample) | PASS; differences from BUSCO's `.faa` only the known +1 residue cases below |
| Rejected hits | 0 in every sample; exported + rejected = hit rows |
| Ortholog counts vs BUSCO `short_summary.json` | equal in every sample |
| Independent verifier (no BuscoPainter, no Biopython code) | 4/4 PASS: expected CDS set = GenBank CDS set; translations re-derived; per-contig files equal the multi-record file |
| Artemis 18.2.0 file reader vs independent parser, every per-contig file | identical locations, strands, joins and `codon_start` |
| Artemis 18.2.0 rendering (offscreen) | joins, frames and trimmed models drawn as expected |
| Features changed vs `T7.1-fixed.2` | only the BUSCO-trimmed models (see below); all others identical |

There are more exported genes than Complete BUSCOs because each copy of a Duplicated BUSCO is
exported, and Fragmented hits are exported as well (exported genes = S rows + D rows + F rows).

### Independent validator

`buscopainter-validate --check-miniprot` runs 25 checks per sample. **All four samples pass all 25
checks.** The checks are:

- **GenBank integrity:** the file parses; record count, IDs and order equal the FASTA; every
  sequence is identical.
- **Features:** all features lie within bounds; strands are valid; locus tags are unique; only
  `source`, `gene` and `CDS` features exist; the gene set equals the CDS set; each CDS lies inside
  its gene on the same strand.
- **Traceability:** every feature traces to a BUSCO row and a GFF mRNA; the gene span equals the
  mRNA span (or the span of the rows BUSCO kept, for a BUSCO-trimmed model); the CDS parts equal
  the GFF CDS rows BUSCO reported, in transcription order; the strand agrees with
  the GFF and with `full_table.tsv`; `codon_start` equals the 5′ phase + 1.
- **Translations:** every attached translation is the clean translation of its CDS, and it agrees
  with miniprot's own protein (`##STA`).
- **Reconciliation:** Missing BUSCOs are never annotated; no row is used twice; every hit is either
  exported or listed in `skipped.tsv`; duplicated copies are distinct and complete.
- **Counts:** ortholog counts equal `short_summary.json`; the painter's summary equals the
  independently recomputed counts; exported + rejected = hit rows.

### Translations compared with miniprot

Each attached translation was compared with miniprot's own protein (`##STA`). In every sample,
the only differences are known miniprot/BUSCO quirks:

- **Stop-less models:** for alignments without a stop codon, miniprot's protein (which BUSCO
  copies into `.faa`) includes 1 extra C-terminal residue, from the codon after the aligned
  region. In all miniprot alignments checked, the GFF CDS length equals the CIGAR span exactly,
  so BuscoPainter keeps the GFF coordinates.
- **BUSCO-trimmed models whose first kept row has phase > 0:** BUSCO's `.faa` has 1 extra
  N-terminal residue, the codon that spans the dropped junction. The GenBank translation starts at
  `/codon_start` and omits it.

There were no other differences.

### Visual inspection

Artemis opens only the first record of a multi-record GenBank file, so an earlier claim that the
multi-record files were "inspected in Artemis" covered only their first contig. Since
T7.1-fixed.2, Artemis 18.2.0's own file reader (`artemis.jar`, run headless) has been run on every
per-contig file: each loads the whole contig with the same features as Biopython reads.

### Whole-assembly Artemis export (after T7.1-fixed.3)

`buscopainter-artemis` (indexed FASTA + tabix-indexed GFF3) was run on the four validated
T7.1-fixed.3 GenBank files. It does not change any gene model.

| Check | Result |
|---|---|
| `buscopainter-validate-artemis` (FASTA copy, `.fai` vs `samtools faidx`, BGZF, GFF3 structure, models vs GenBank and `features.tsv`, phases, per-contig counts, tabix queries) | 4/4 PASS; 19,091 gene models, 1,849 multi-segment |
| Artemis 18.2.0 GUI, scripted (`tools/artemis_gui_check/`): command-line route, every contig chosen in the drop-down list, feature view scrolled along each contig | 4/4 PASS: 725/725 contigs offered with the right length; every model loaded with the same coordinates, strand and joins as the GenBank file, and drawn |
| Same, menu route (File → Open, File → Read An Entry), one assembly | PASS |
| Multi-record GenBank in the Artemis GUI | first contig only (as documented) |
| Uncompressed GFF3, or BGZF GFF3 without `.tbi` | not usable in Artemis (docs/artemis.md) |

### Negative control

The pre-fix script (T7.1) **exits with code 0 on the validation data but fails 7 validator
checks**. It wrote more genes than there were hits, shifted CDS boundaries by the GFF phase,
produced non-unique locus tags, and reported the Duplicated hit count as if it were the ortholog
count. This is why exit code 0 is not taken as evidence of correctness.

## T7.1-fixed.3 change (BUSCO-trimmed models)

In under 1 % of models, BUSCO's `full_table.tsv` span covers only some whole CDS rows of the
miniprot model, and BUSCO's own `.faa` protein is the translation of those rows only. Earlier
versions (T5, T7.1, T7.1-fixed.1/.2) painted the whole GFF model, including segments BUSCO had
dropped, up to 10.8 kb from the gene. T7.1-fixed.3 paints BUSCO's rows only and lists the dropped
rows in a `/note`. All other features are identical to T7.1-fixed.2.

## T7.1-fixed.2 change (reporting and Artemis)

T7.1-fixed.2 added per-contig export, a scope `/note` on every `source` feature, rejection of
partial GFF models and exit code 3 on a failed reconciliation. On the validation data, every gene
and CDS feature and every TSV report was identical to T7.1-fixed.1. The GenBank text differed only
by the new `/note`.

## Packaging equivalence

The installed `buscopainter` command was run on the four samples and compared with the validated
reference outputs: GenBank and TSV outputs were byte-identical (Python 3.14 and a fresh Python 3.11
environment, checked at T7.1-fixed.1), and `--dry-run` wrote the same reports and no GenBank.

## Automated tests

```bash
pip install -r requirements-lock.txt -e ".[test]"
python -m pytest
```

There were 45 tests at T7.1-fixed.1 (79 now, including the 12 T7.1-fixed.2 and 8 T7.1-fixed.3
regression tests and 14 tests of the Artemis export, 7 of which need `tabix` and are skipped
without it), all on synthetic data with known answers. They pass on Python 3.11, 3.12, 3.13
and 3.14 (installed from the built wheel into fresh environments). They cover:

- plus- and minus-strand coordinates;
- GFF phase 1 with `/codon_start=2` and no boundary shift;
- multi-segment `complement(join())` CDS in transcription order;
- three Duplicated copies (two on one contig) kept distinct;
- Fragmented hits without a stop codon;
- Missing rows, including 2-field rows that come first in the table;
- no translation for frame-inconsistent, internal-stop or miniprot-frameshift models;
- every rejection type, with no clamping;
- models that end exactly at a contig end (accepted) and 1 bp past it (rejected);
- ortholog versus hit-row counts;
- unique locus tags; one gene + one CDS per hit;
- `--dry-run` writes the same reports and no GenBank;
- consistency between the GenBank, features, summary and issues files;
- the synthetic example reproduces `examples/expected_output/` byte for byte and passes the
  validator;
- the CLI entry points;
- the SHA-256 of the reference implementation;
- the Artemis export: designed coordinates and phases, escaping, sorting, a refused mismatched
  FASTA, and the validator failing on six deliberate corruptions (coordinate, strand, phase,
  missing model, changed FASTA base, missing `.tbi`).

## Reproducing the validation on your own data

```bash
buscopainter --config config.json --outdir results
buscopainter-validate --config config.json --outdir results --check-miniprot --report results/validation.json
```

`--check-miniprot` requires `busco_dir` to be the parent of `run_*` (where BUSCO's `logs/` folder
is).
