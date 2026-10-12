# Input formats

## What has been tested

| Item | Tested | Not tested |
|---|---|---|
| BUSCO | 6.0.0, genome mode | BUSCO 5.x, protein and transcriptome modes |
| Gene predictor | miniprot 0.18-r281 (BUSCO's genome-mode default) | MetaEuk (`--metaeuk`), AUGUSTUS (`--augustus`) |
| Lineage | `euglenozoa_odb12` (n = 3,611) | other lineages, auto-lineage runs |
| Assemblies | hifiasm primary assemblies, 44–333 contigs | chromosome-scale assemblies with very long IDs, highly fragmented short-read assemblies |

The code reads the general BUSCO v5/v6 layout, which is the same for all lineages, so it may work
outside the tested range. Run `--dry-run` and the validator to check.

## Assembly FASTA

- It must be **the same assembly BUSCO was run on**. BuscoPainter cannot detect a different file
  with the same contig names. To check, compare its checksum with the `in` path in BUSCO's
  `short_summary.json` (`parameters.in`).
- The sequence ID is the first word of each header line (Biopython convention). IDs must be unique
  and must equal the `Sequence` column of `full_table.tsv`.
- Plain-text FASTA only. Compressed FASTA has not been tested.
- Every FASTA record appears in the GenBank output, including contigs without BUSCO hits.

## BUSCO run directory

```text
<busco_dir>/                       # pass this or run_<lineage>/ as busco_dir
├── logs/miniprot_align_<lineage>_out.log     # only for buscopainter-validate --check-miniprot
└── run_<lineage>/
    ├── full_table.tsv                         # or full_table.tsv.gz (painter only)
    ├── short_summary.json                     # optional; used to cross-check counts
    └── busco_sequences/
        ├── single_copy_busco_sequences/<id>.gff
        ├── multi_copy_busco_sequences/<id>.gff
        └── fragmented_busco_sequences/<id>.gff
```

If `busco_dir` contains `full_table.tsv`, it is used as the run directory. Otherwise it must contain
exactly one `run_*` subdirectory with a `full_table.tsv`. The lineage label in the notes is taken
from the `run_<lineage>` name.

## `full_table.tsv`

BUSCO v5/v6 format. Lines beginning with `#` are skipped. The columns are:

```text
Busco id  Status  Sequence  Gene Start  Gene End  Strand  Score  Length  OrthoDB url  Description
```

- `Status` must be `Complete`, `Duplicated`, `Fragmented` or `Missing`. Any other value is a
  fatal error.
- `Missing` rows can have just two fields. They are counted but never placed on a contig.
- Every non-Missing row needs `Sequence`, integer `Gene Start` and `Gene End`, and a `+` or `-`
  strand.
- A Duplicated BUSCO has one row per copy.

**Coordinate convention.** In the tested BUSCO 6 / miniprot data, `Gene Start` and `Gene End` are
written in strand order (the higher number first on the minus strand). They use a 0-based lower
bound and exclude the stop codon. The miniprot GFF uses 1-based inclusive coordinates and includes
the stop codon. **BuscoPainter takes all coordinates from the GFF.** It uses `full_table.tsv` only
to decide which GFF model belongs to which hit.

## Per-BUSCO GFF (miniprot)

BUSCO writes one GFF per BUSCO id. This file holds **every** copy of that BUSCO. BuscoPainter uses
these rows:

| Type | Required attributes | Use |
|---|---|---|
| `mRNA` | `ID` (unique within the file) | Gene span; `Identity`, `Target`, `Frameshift`, `StopCodon` are read if present |
| `CDS` | `Parent` = an mRNA `ID`; phase 0, 1 or 2 | CDS segments |
| `stop_codon` | `Parent` | Checked to lie inside the mRNA span |

Other feature types are ignored. A model is **rejected** (`INVALID_GFF_MODEL`) if any of these is
true:

- a child row is on a different contig or strand from its mRNA, or lies outside the mRNA span;
- two CDS rows overlap;
- the mRNA has no CDS rows.

Malformed lines (wrong column count, non-integer coordinates, invalid strand or phase, orphan
children) are reported as `GFF_PARSE_PROBLEM` warnings.

## How hits are matched to GFF models

For each BUSCO id, every non-Missing `full_table.tsv` row is paired with one GFF `mRNA`. The pair
must be on the same contig with overlapping spans. The closest span wins, measured as
|Δstart| + |Δend|. Each mRNA is used at most once. When two equally close models are free, the hit
is rejected as `AMBIGUOUS_MATCH` rather than guessed. GFF models that no row claims are reported as
`INFO GFF_MODEL_NOT_IN_FULL_TABLE` and are not exported.

If the `full_table.tsv` span differs from the GFF mRNA span by more than 3 bp at either end, the
GFF span is exported, the hit is flagged `FULLTABLE_SPAN_DIFFERS_FROM_GFF`, and no translation is
attached.

The GFF strand must equal the `full_table.tsv` strand. Otherwise the hit is rejected as
`STRAND_MISMATCH`.
