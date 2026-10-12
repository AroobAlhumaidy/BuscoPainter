---
name: Bug report
about: Wrong output, a crash, or a hit that was rejected unexpectedly
labels: bug
---

**What happened**
<!-- One or two sentences. -->

**Versions**
- BuscoPainter: <!-- python -c "import buscopainter; print(buscopainter.__version__)" -->
- Python:
- BUSCO version and gene predictor (miniprot / metaeuk / augustus):
- Lineage dataset (e.g. eukaryota_odb12):

**Command and config**
```
buscopainter --config ... --outdir ...
```

**Evidence**
- The relevant lines of `<prefix>.busco_painter.log` (look for `WARNING` and `ERROR`).
- The relevant rows of `<prefix>_<sample>.annotation_issues.tsv`.
- If possible, the `buscopainter-validate` output.
- For a single wrong feature: the `full_table.tsv` row and the matching GFF lines.

Please do not attach unpublished genomes. A few lines of the files above are usually enough.
If they are not, a small synthetic reproduction is best.
