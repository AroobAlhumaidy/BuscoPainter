# Output formats

All files are written to `--outdir`. `<prefix>` comes from `--prefix` or the config, and `<sample>`
is the sample `name`.

## GenBank: `<prefix>_<sample>.gbk`

The file has one record per FASTA sequence, in FASTA order. The sequence is copied unchanged. The
`LOCUS` date is Biopython's default (`01-JAN-1980`) and the division is `UNK`.

| Feature | Location | Qualifiers |
|---|---|---|
| `source` | whole record | `/organism`, `/mol_type="genomic DNA"`, `/isolate` (if given), `/note` stating that only BUSCO `<lineage>` marker genes are annotated |
| `gene` | GFF `mRNA` span (or, for a BUSCO-trimmed model, the span of the CDS rows BUSCO kept) | `/locus_tag`, `/note` (BUSCO id, status, GFF evidence) |
| `CDS` | GFF `CDS` rows; one segment, or `join(...)` in transcription order; `complement(...)` on the minus strand | `/locus_tag`, `/product`, `/codon_start`, `/transl_table=1`, several `/note`, and `/translation` only when it is clean |

There is exactly one `gene` and one `CDS` per exported BUSCO hit, and they share the same
`locus_tag`. The tag has the form `<sample>_<busco id>_<miniprot mRNA ID>`, and it is unique
within a sample.

### Coordinates and phase

GFF coordinates (1-based, inclusive) are written unchanged. **GFF phase never moves a boundary.**
The phase of the 5′-most CDS segment becomes `/codon_start` (phase + 1). For example, a CDS with
phase 1 gets `/codon_start=2`, and its translation starts at the second base.

### CDS notes

| Note | When |
|---|---|
| `BUSCO <id> (<status>), lineage <lineage>; full_table line <n>: <contig>:<start>-<end>(<strand>) score=… length=…` | always |
| `GFF evidence: <gff file> mRNA <ID> <contig>:<start>-<end>(<strand>); Identity=…; Target=…` | always |
| `CDS joins <k> miniprot segments separated by [gaps] bp; cis-introns are rare in trypanosomatids, so this structure is unverified` | multi-segment models. The wording is fixed and mentions trypanosomatids whatever the organism. |
| `GFF phases are not consistent with one reading frame across segments (possible frameshift); exported as given` | frame-inconsistent models |
| `translation omitted: <reason>` | no `/translation` attached |

### When `/translation` is attached

All of the following must hold:

1. The GFF phases of all segments form one continuous reading frame.
2. BUSCO's `full_table.tsv` span covers the exported CDS (within 3 bp at each end). A model that
   BUSCO trimmed to whole CDS rows (`FULLTABLE_SPAN_TRIMMED`) meets this after trimming.
3. The miniprot `mRNA` and its `CDS` rows have no `Frameshift` or `StopCodon` attribute.
4. The CDS length after `codon_start` is a multiple of 3.
5. Translating with the standard code (table 1) gives no internal stop codon, and the protein is
   not empty or made only of `X`.

The translation is a literal translation of the CDS. A terminal stop is removed. **The first codon
is not converted to M**, so an alternative start codon such as GTG appears as V.

## `<prefix>_<sample>.features.tsv`

One row per exported hit.

| Column | Meaning |
|---|---|
| `sample`, `locus_tag`, `busco_id`, `status`, `contig`, `contig_length`, `strand`, `mrna_id` | Identity of the hit |
| `gene_start`, `gene_end` | gene span written to the GenBank (GFF mRNA span, or the kept CDS rows of a BUSCO-trimmed model), 1-based inclusive |
| `cds_location` | GenBank-style location, 1-based |
| `cds_segments` | Number of GFF CDS rows joined |
| `codon_start` | 5′ segment phase + 1 |
| `frame_consistent` | `True` if the phases form one reading frame |
| `translation_attached`, `translation_length`, `translation_issue` | Translation outcome and the reason if it was omitted |
| `starts_with_ATG` | `True` only if `codon_start` is 1 and the CDS starts with ATG |
| `ends_with_stop` | `True` if the translated CDS ends in a stop codon. **It is only evaluated when translation was attempted.** For models blocked earlier (frameshift, partial BUSCO span, frame inconsistency, length not a multiple of 3), it is `False` without having been checked. |
| `flags` | `;`-separated: `MULTI_SEGMENT_MODEL`, `FRAME_INCONSISTENT`, `FULLTABLE_SPAN_TRIMMED`, `FULLTABLE_SPAN_DIFFERS_FROM_GFF`, `MINIPROT_FRAMESHIFT`, `MINIPROT_STOPCODON`, `NO_TRANSLATION` |
| `ft_line`, `ft_start`, `ft_end`, `ft_score`, `ft_length` | The `full_table.tsv` line number and values, as given (strand order, BUSCO convention) |
| `gff_file` | GFF path relative to the run directory |
| `product` | `/product` value |

## `<prefix>_<sample>.annotation_issues.tsv` and `.skipped.tsv`

Columns: `sample, busco_id, status, sequence, ft_start, ft_end, ft_strand, mrna_id, severity, code, detail`.
`skipped.tsv` contains only the `REJECTED` rows. If it holds only a header line, nothing was
rejected.

| Severity | Code | Meaning |
|---|---|---|
| REJECTED | `MISSING_GFF` | No `<id>.gff` in any status folder |
| REJECTED | `NO_OVERLAPPING_GFF_MODEL` | No GFF mRNA on that contig overlaps the hit |
| REJECTED | `MODEL_ALREADY_ASSIGNED` | All overlapping models were claimed by other rows |
| REJECTED | `AMBIGUOUS_MATCH` | Two or more equally close models; not guessed |
| REJECTED | `INVALID_GFF_MODEL` | Inconsistent GFF model (see input_formats.md) |
| REJECTED | `STRAND_MISMATCH` | GFF strand differs from `full_table.tsv` |
| REJECTED | `MISSING_CONTIG_IN_FASTA` | Contig not in the FASTA |
| REJECTED | `OUT_OF_RANGE` | Model extends past the contig end (never clamped) |
| REJECTED | `CONFLICTING_STATUS` | One BUSCO id with more than one status |
| REJECTED | `DUPLICATE_LOCUS_TAG` | Locus tag collision |
| WARNING | `MULTI_SEGMENT_MODEL` | ≥ 2 CDS segments; the gap sizes are given |
| WARNING | `FRAME_INCONSISTENT` | Phases are not one reading frame |
| WARNING | `FULLTABLE_SPAN_TRIMMED` | BUSCO's span is a run of whole CDS rows of a longer miniprot model (BUSCO dropped the others). Only BUSCO's rows are exported; the dropped rows are listed |
| WARNING | `FULLTABLE_SPAN_DIFFERS_FROM_GFF` | Spans differ by > 3 bp and do not fall on whole CDS rows; the GFF span is exported without a translation |
| WARNING | `TRANSLATION_OMITTED` | Reason in `detail` |
| WARNING | `MULTIPLE_ROWS_FOR_SINGLE_COPY` | A Complete BUSCO has more than one row |
| WARNING | `GFF_IN_UNEXPECTED_FOLDER` | GFF found under a different status folder |
| WARNING | `GFF_PARSE_PROBLEM` | Malformed GFF line |
| INFO | `GFF_MODEL_NOT_IN_FULL_TABLE` | GFF model not claimed by any row; not exported |

A model can have several warnings, so the number of warning events is usually larger than the
number of flagged models.

## `<prefix>_<sample>.contig_summary.tsv`

One row per contig plus a `TOTAL` row: `contig, length, Complete_hit_rows, Duplicated_hit_rows,
Fragmented_hit_rows, exported_genes, exported_CDS, cds_segments`. Missing BUSCOs have no location,
so they do not appear here.

## `<prefix>.busco_summary.tsv`

One row per sample.

| Column | Meaning |
|---|---|
| `busco_{S,D,F,M}_orthologs`, `busco_C_orthologs`, `busco_n_orthologs` | **Distinct BUSCO ids** per status: BUSCO's own convention, C = S + D |
| `busco_{S,D,F}_hit_rows` | **Rows** in `full_table.tsv` per status. A Duplicated BUSCO has ≥ 2 rows, one per copy. |
| `short_summary_agrees` | Ortholog counts equal `short_summary.json` (`True`/`False`), or empty if that file is absent |
| `exported_genes`, `exported_CDS` | One per exported hit |
| `cds_segments`, `multi_segment_models`, `translations_attached` | Totals over exported hits |
| `rejected_hits`, `warnings` | Number of REJECTED hits and WARNING events |
| `dry_run` | Whether GenBank was skipped |
| `reconciled` | `exported_genes + rejected_hits` equals the S + D + F hit rows (`False` → exit code 3) |

`exported_genes + rejected_hits` equals the sum of the S, D and F hit rows.

## Duplicated orthologs versus duplicated hits

BUSCO reports D as a number of **orthologs** (BUSCO ids found in ≥ 2 complete copies). The
`full_table.tsv` has one row per **copy**. For example, a BUSCO with three copies counts once in
D but has three Duplicated rows. BuscoPainter exports every copy as its own feature, and it reports both numbers
in separate, labelled columns.

## Per-contig GenBank (optional)

With `--per-contig` or `export_per_contig: true`, the painter also writes
`<prefix>_<sample>_per_contig/<prefix>_<sample>_contig-<id>.gbk`, one single-record file per contig.
Characters other than letters, digits, `.`, `_` and `-` in a contig id become `_` in the file name;
if two ids map to the same name the run stops with an error. Records are identical to those in the
multi-record file. Use these files in Artemis, which opens only the first record of a multi-record
GenBank file. Validated on every contig of the four validation assemblies (docs/validation.md).
