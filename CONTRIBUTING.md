# Contributing

## Reporting bugs and requesting features

Please open a GitHub issue: <https://github.com/AroobAlhumaidy/BuscoPainter/issues>.

For a bug, include:

- the BuscoPainter version (`python -c "import buscopainter; print(buscopainter.__version__)"`),
  the Python version, the BUSCO version, the gene predictor and the lineage;
- the command and the config (you can redact the paths);
- the relevant `WARNING` and `ERROR` lines from the log and rows from `annotation_issues.tsv`;
- the `buscopainter-validate` output, if you ran it.

**Please do not attach unpublished genomes or full BUSCO runs.** A few lines of `full_table.tsv`
and the GFF are usually enough. A small synthetic reproduction (see `examples/make_example_data.py`)
is ideal.

For a feature request, describe the use case and the input data type (BUSCO version, mode, gene
predictor).

## Changing the code

`src/buscopainter/painter.py` is the validated scientific reference. Any change to it can change
the output, so a change must:

1. come with a test on synthetic data with a known answer;
2. keep `python -m pytest` passing. `test_reference_implementation_unchanged` fails by design
   until its SHA-256 is updated on purpose;
3. be listed in `CHANGELOG.md`, saying whether it changes the output;
4. if the output changes, be validated again with `buscopainter-validate` on real data, with a
   summary of the result recorded in `docs/validation.md`. Do not commit the real data or
   per-sample results from unpublished genomes.

`src/buscopainter/validate.py` must stay independent. It must not import the painter.

## Versioning

Versions follow Semantic Versioning (`MAJOR.MINOR.PATCH`). A change in the scientific output is at
least a MINOR bump while the version is below 1.0, and it is always listed in the changelog.

## Licence

BuscoPainter is released under the MIT License (see `LICENSE`). By contributing, you agree that
your contribution is released under the same licence.
