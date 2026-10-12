# Synthetic example

Everything in `data/` is **invented**: random DNA with hand-placed open reading frames, plus
BUSCO-style files (`full_table.tsv`, `short_summary.json`, miniprot-style GFFs) made by
`make_example_data.py`. It is not real BUSCO output and has no biological meaning. It is safe to
redistribute.

## Run it

Run these commands from this directory, because the config uses relative paths:

```bash
cd examples
buscopainter --config config.json --outdir my_output
buscopainter-validate --config config.json --outdir my_output
```

The outputs in `my_output/` should be byte-identical to `expected_output/`, except for the log
file, which holds timestamps.

To regenerate `data/` (the result is deterministic):

```bash
python make_example_data.py
```

## What it contains

There are 3 contigs (`contig_3` has no hits) and 9 BUSCO ids: 5 Complete, 1 Duplicated (3 copies),
1 Fragmented and 2 Missing. That gives 9 hit rows, all exported.

| BUSCO | Case | Expected result |
|---|---|---|
| `ex1at0` | + strand, one segment, with stop | CDS `201..386`, translation attached |
| `ex2at0` | − strand, two segments, frame-consistent | `complement(join(801..996,1037..1086))`, translation attached, `MULTI_SEGMENT_MODEL` |
| `ex3at0` | GFF phase 1 | CDS `1501..1627` (not shifted), `/codon_start=2`, translation attached |
| `ex4at0` | Duplicated: 2 copies on `contig_1`, 1 on `contig_2` (− strand) | 3 separate features with distinct locus tags |
| `ex5at0` | Fragmented, no stop codon | exported, translation attached |
| `ex6at0` | miniprot `Frameshift=1` | exported, `translation omitted` |
| `ex7at0` | two segments with inconsistent phases | `FRAME_INCONSISTENT`, `MULTI_SEGMENT_MODEL`, `translation omitted` |
| `ex8at0`, `ex9at0` | Missing (`ex9at0` is the first row of the table) | counted (M = 2), never placed |

Summary (`expected_output/example.busco_summary.tsv`): S = 5, D = 1 (3 hit rows), F = 1, M = 2,
C = 6, n = 9. The run exports 9 genes and 9 CDS (11 segments) and attaches 7 translations. It
rejects nothing and logs 5 warning events on 3 models.
