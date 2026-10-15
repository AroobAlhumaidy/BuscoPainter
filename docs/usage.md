# Usage

## Commands

| Command | What it runs |
|---|---|
| `buscopainter` (or `python -m buscopainter`) | The painter (`buscopainter.painter`) |
| `buscopainter-validate` (or `python -m buscopainter.validate`) | The independent output validator (`buscopainter.validate`) |
| `buscopainter-artemis` (or `python -m buscopainter.artemis_export`) | Optional: indexed FASTA + GFF3 so Artemis can browse every contig ([artemis.md](artemis.md)) |
| `buscopainter-validate-artemis` (or `python -m buscopainter.validate_artemis`) | Checks the `buscopainter-artemis` output |

## `buscopainter` options

| Option | Default | Meaning |
|---|---|---|
| `--config FILE` | required | JSON config, or YAML if PyYAML is installed (`.yaml` / `.yml`) |
| `--outdir DIR` | config `outdir`, else `results` | Output directory, created if needed |
| `--prefix TEXT` | config `prefix`, else `Project` | Prefix for every output file name |
| `--samples NAME ...` | all | Process only these sample names |
| `--dry-run` | off | Parse, match and check everything, and write the TSV reports, but write no GenBank |
| `--per-contig` | off | Also write one single-record GenBank per contig (for Artemis); same as `export_per_contig` |
| `--log-level LEVEL` | `INFO` | Log level for the log file (`DEBUG`, `INFO`, `WARNING`, ...) |

The `--help` text shows the internal version string `T7.1-fixed.3`. This is the validated reference
implementation that release 0.1.0 packages (gene models identical to `T7.1-fixed.1` except for
BUSCO-trimmed models; see CHANGELOG).

## Config file

```json
{
  "prefix": "myproject",
  "organism": "Genus species",
  "topology": "linear",
  "export_per_contig": false,
  "use_orthodb_api": false,
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

Top-level keys:

| Key | Default | Meaning |
|---|---|---|
| `samples` | required | List of samples (below) |
| `prefix` | `Project` | Output file prefix (`--prefix` overrides it) |
| `outdir` | `results` | Output directory (`--outdir` overrides it) |
| `organism` | `Unknown` | Written to the GenBank `ORGANISM` line and the `/organism` source qualifier |
| `topology` | `linear` | GenBank topology for every record (`linear` or `circular`) |
| `export_per_contig` | `false` | Also write one GenBank file per contig, in `<prefix>_<sample>_per_contig/` (use these in Artemis) |
| `use_orthodb_api` | `false` | Query the OrthoDB web API for product names when the full_table `Description` is empty. Needs `requests`. **Not used in the validation.** |

Per-sample keys:

| Key | Required | Meaning |
|---|---|---|
| `fasta` | yes | Assembly FASTA (the file BUSCO was run on) |
| `busco_dir` | yes | The BUSCO `run_<lineage>/` directory, or its parent if that holds exactly one `run_*` directory with a `full_table.tsv` |
| `name` | recommended | Sample name, used in file names and locus tags. Defaults to the FASTA file stem. |
| `isolate` | no | `/isolate` on the source feature |
| `lineage` | no | Lineage label for the notes. By default it is taken from the `run_<lineage>` folder name. |
| `organism`, `topology`, `prefix`, `export_per_contig` | no | Per-sample overrides |

Rules:

- **Relative paths are resolved from the current working directory**, not from the config file.
- The config path, the FASTA paths and the BUSCO paths **must not contain whitespace**. If they do,
  the run stops with exit code 2.
- `buscopainter-validate` reads **JSON configs only**.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | All samples were processed. **This does not mean every hit was exported or is correct.** |
| 2 | Fatal input error, for example a missing file, an unresolvable `run_*` directory, duplicate FASTA IDs, a malformed `full_table.tsv`, or whitespace in a path. Paths, BUSCO directories and full tables for all samples are checked before any sample is processed. FASTA problems (empty file, duplicate IDs) are found only when that sample is processed, so files for earlier samples may already be written. |
| 3 | Count reconciliation failed: `exported + rejected` does not equal the number of non-Missing hit rows (see below) |
| other | Unexpected error (Python traceback) |

After a run:

1. Check the log for `ERROR` lines. If `exported + rejected` does not equal the number of non-Missing
   hit rows, the painter logs an `ERROR`, sets `reconciled` to `False` in the summary and exits
   with code 3.
2. Check `<prefix>_<sample>.skipped.tsv` for rejected hits.
3. Check `<prefix>.busco_summary.tsv`. If `short_summary_agrees` is `False`, the `full_table.tsv`
   counts disagree with BUSCO's `short_summary.json`.
4. Run the validator.

## Recommended workflow

```bash
buscopainter --config config.json --outdir results --dry-run
# read results/*.annotation_issues.tsv and results/*.skipped.tsv
buscopainter --config config.json --outdir results
buscopainter-validate --config config.json --outdir results --report results/validation.json
```

## `buscopainter-validate` options

| Option | Meaning |
|---|---|
| `--config FILE` | Same JSON config as the painter |
| `--outdir DIR` | Painter output directory to check |
| `--prefix TEXT` | Must match the painter's prefix (default: config `prefix`, else `Project`) |
| `--samples NAME ...` | Check only these samples |
| `--check-miniprot` | Also compare each attached translation with miniprot's own protein (`##STA` lines in `<busco_dir>/logs/miniprot_align_*_out.log`). Here `busco_dir` must be the parent of `run_*`. |
| `--check-faa` | Also compare translations with BUSCO's `busco_sequences/*/<id>.faa` |
| `--report FILE` | Write a JSON report |

It prints `OVERALL: PASS` or `OVERALL: FAIL`, and exits with 0 only if every check passes for
every sample (otherwise 1). The validator does not import the
painter. It re-implements parsing, matching and counting separately.

The validator expects the standard layout (`run_*/full_table.tsv`, uncompressed, and
`run_*/busco_sequences/<status folder>/<id>.gff`). The painter also accepts `full_table.tsv.gz`
and some older folder layouts, but the validator does not.

## `buscopainter-artemis` options

Reads the painter's GenBank output; it does not read BUSCO and does not change any gene model.
See [artemis.md](artemis.md) for how to open the result.

| Option | Meaning |
|---|---|
| `--config FILE` | Painter config (JSON). Exports every sample (or `--samples`) |
| `--outdir DIR` | With `--config`: the painter output directory. With `--gbk`: the export directory |
| `--export-root DIR` | With `--config`: write the `<prefix>_<sample>_artemis/` folders here instead of `--outdir` |
| `--prefix TEXT` | Painter prefix (default: config `prefix`, else `Project`) |
| `--samples NAME ...` | Only these samples |
| `--gbk FILE --fasta FILE` | Convert one painter GenBank file and the FASTA it was made from |
| `--tabix PATH` | `tabix` program (default: search `PATH`). Without it the command stops with exit code 2, because Artemis needs the `.tbi` index |
| `--report FILE` | Write a JSON summary |

Exit codes: 0 done; 2 error (for example, `tabix` is missing, the FASTA does not hold exactly the
GenBank sequences, or a FASTA record has uneven line lengths so it cannot be indexed).

## `buscopainter-validate-artemis` options

`--config`, `--outdir`, `--export-root`, `--prefix`, `--samples`, `--tabix` and `--report`, as
above. It prints `OVERALL: PASS` or `OVERALL: FAIL` and exits with 0 only on PASS.
