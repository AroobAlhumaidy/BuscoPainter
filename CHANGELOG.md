# Changelog

Versions follow [Semantic Versioning](https://semver.org/). Any change that can alter the
scientific output (coordinates, features, translations, flags or counts) is listed explicitly,
and the release is validated again before it is published.

## [0.1.0] - unreleased

The first public release. It packages internal version `T7.1-fixed.3` as `buscopainter.painter`:
the validated reference implementation `buscoPainter.T7.1.fixed.py` (`T7.1-fixed.1`) plus the
changes under "T7.1-fixed.2" and "T7.1-fixed.3" below.

### Release
- **Licence** MIT (`LICENSE`); set in `pyproject.toml` and `CITATION.cff`.
- **Repository** renamed to `BuscoPainter`; URLs updated in the README, `CONTRIBUTING.md`,
  `pyproject.toml` and `CITATION.cff`.
- **Docs** the README states that BuscoPainter transfers BUSCO marker annotations and is not a
  gene prediction tool. Per-sample results from the unpublished validation assemblies were removed
  from the public documentation; `docs/validation.md` gives a summary.

### T7.1-fixed.3 (2026-10-12), final verification
**Changes scientific output for BUSCO-trimmed models** (under 1 % of models in the validation
data); every other feature is identical to `T7.1-fixed.2`.
- **Fixed** when BUSCO 6 drops terminal miniprot segments, `full_table.tsv` and BUSCO's `.faa`
  protein cover only the remaining whole CDS rows, but the per-BUSCO GFF keeps the whole model.
  Earlier versions (T5, T7.1, T7.1-fixed.1/.2) painted the whole GFF model, including segments
  BUSCO had dropped, up to 10.8 kb from the gene. Now only BUSCO's rows are exported, the gene
  span is the span of those rows, a `/note` lists the dropped rows, and the event is logged as
  `FULLTABLE_SPAN_TRIMMED`. A span that does not fall on whole rows keeps the old behaviour
  (`FULLTABLE_SPAN_DIFFERS_FROM_GFF`, no translation). Most trimmed models now get a
  translation (the others keep a stated reason: frameshift or length not a multiple of 3), and
  each agrees with BUSCO's `.faa` protein.
- **Validator** expects BUSCO's rows for trimmed models, and accepts BUSCO's extra N-terminal
  junction residue (`.faa`) and miniprot's untrimmed `##STA` protein for them.
- **Tests** 8 regression tests (`tests/test_regressions_t7_1_fixed_3.py`): trimming on both
  strands, a kept segment starting with phase 2, a frameshift only in a dropped segment, the
  no-trim fallback, and the validator accepting trimmed and rejecting untrimmed output.

### T7.1-fixed.2 (2026-10-12), after a review of the GenBank output in Artemis
No gene model, coordinate, translation or count changes on the validation data (every gene and
CDS feature identical to `T7.1-fixed.1`).
- **Added** `--per-contig` CLI flag (same as `export_per_contig`). Artemis opens only the first
  record of a multi-record GenBank file; the per-contig files are the ones to open in Artemis.
  They are now written to `<prefix>_<sample>_per_contig/` (previously, untested, into the output
  directory itself), with file-name-safe contig ids and a collision check.
- **Added** a `/note` on every record's `source` feature: only BUSCO `<lineage>` marker genes are
  annotated; this is not a whole-genome annotation.
- **Added** a log line, when per-contig export is off, saying Artemis shows only the first record.
- **Fixed** a GFF CDS row that failed parsing (bad phase, interval or strand) was dropped while its
  model was still exported with a segment missing (warning only). The model is now rejected as
  `INVALID_GFF_MODEL`. No such rows exist in the validation data.
- **Changed** a failed count reconciliation now gives exit code 3 (was 0 with an `ERROR` line) and
  a `reconciled` column in `busco_summary.tsv`.
- **Docs** corrected the claim that the four multi-record GenBank files were "inspected in
  Artemis"; documented lineage scope (non-marker ortholog groups, repeat proteins).
- **Tests** 12 regression tests (`tests/test_regressions_t7_1_fixed_2.py`): 3-segment models with
  every internal phase on both strands, duplicated copies on opposite strands, GFF models absent
  from `full_table.tsv`, partial models, scope note, per-contig export, exit code 3.

### Added
- `buscopainter` and `buscopainter-validate` commands, and `python -m buscopainter`.
- A synthetic example with expected outputs, and documentation in `docs/`.
- 21 release tests and 20 regression tests (65 in total): packaging, CLI, example reproducibility, translation integrity,
  contig boundaries, output consistency and dry-run equivalence.

### Fixed relative to the earlier internal script T7.1
- The GFF phase no longer moves CDS boundaries. It sets `/codon_start` instead.
- There is one `gene` and one `CDS` per BUSCO hit. A multi-segment model becomes one `join()` CDS.
- Duplicated copies on the same contig are matched one-to-one and are no longer cross-assigned.
- Locus tags are unique (`<sample>_<buscoid>_<mRNA ID>`).
- Out-of-range coordinates are rejected and reported, never clamped.
- Strand mismatches and invalid GFF models are rejected and reported.
- A `/translation` is attached only when it is clean. Otherwise the reason is noted.
- Ortholog and hit-row counts are reported in separate, labelled columns.
- `--dry-run` now writes no GenBank (before, it wrote one).
- `full_table.tsv` parsing no longer fails when a 2-field Missing row comes first.
- Every rejection and warning is written to `annotation_issues.tsv`.
