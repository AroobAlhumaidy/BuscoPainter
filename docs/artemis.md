# Viewing a whole assembly in Artemis

Artemis opens **only the first record** of a GenBank file that holds several. So the main
`<prefix>_<sample>.gbk` shows one contig in Artemis. That file is still correct GenBank for other
tools.

To browse **every contig** in one Artemis window, use the indexed export. It writes:

- the assembly FASTA (a byte-identical copy) with a samtools `.fai` index, and
- one GFF3 file holding the BUSCO genes of all contigs, BGZF-compressed, with a tabix `.tbi`
  index.

Artemis then shows a **contig drop-down list** above the sequence view. Choose a contig and its
BUSCO genes appear. The sequences are not joined together or changed in any way.

Tested with Artemis 18.2.0 (macOS) on four assemblies (725 contigs in total): every contig could
be chosen and showed exactly the genes in the GenBank file (see [Validation](#validation)).

## What you need

- Artemis 18 or later.
- `tabix` from [htslib](https://www.htslib.org/) (for example `conda install -c bioconda htslib`
  or `samtools`). **Without the `.tbi` index this view does not work** (see
  [What does not work](#what-does-not-work)).

You do not need `samtools` or `bgzip`. BuscoPainter writes the `.fai` index and the BGZF file
itself.

## Commands

```bash
# 1. Run the painter as usual (writes results/myproject_S1.gbk and friends)
buscopainter --config config.json --outdir results

# 2. Write the Artemis files for every sample
buscopainter-artemis --config config.json --outdir results

# 3. Optional: check them
buscopainter-validate-artemis --config config.json --outdir results
```

Step 2 makes one folder per sample, for example `results/myproject_S1_artemis/`:

| File | What it is |
|---|---|
| `S1.fasta` | Exact copy of the assembly FASTA named in the config |
| `S1.fasta.fai` | FASTA index (same as `samtools faidx` writes) |
| `myproject_S1.busco.gff3.gz` | BUSCO genes on all contigs, BGZF-compressed GFF3 |
| `myproject_S1.busco.gff3.gz.tbi` | tabix index |
| `README.txt` | The opening steps below |

Options:

- `--export-root DIR` writes the `_artemis` folders to `DIR` instead of the painter folder.
- `--gbk FILE --fasta FILE --outdir DIR` converts one GenBank file (no config needed).
- `--tabix PATH` if `tabix` is not on your `PATH`.

If `tabix` is not found, `buscopainter-artemis` stops with exit code 2. Install htslib/tabix, or
pass the executable with `--tabix /path/to/tabix`.

## Open it in Artemis

**From the menus** (macOS app, Windows, Linux):

1. **File -> Open ...** and choose `S1.fasta`.
2. Artemis asks *"Indexed FASTA Found"*. Click **Use index**.
3. In the new window, **File -> Read An Entry ...** and choose `myproject_S1.busco.gff3.gz`.
4. Artemis asks *"Indexed GFF"*. Click **Use index**.
5. Choose a contig in the drop-down list next to the FASTA name, at the top of the window.

**From a terminal**, if the `art` launcher is installed:

```bash
cd results/myproject_S1_artemis
art S1.fasta + myproject_S1.busco.gff3.gz
```

On macOS the app does not install `art`. This does the same:

```bash
A=/Applications/Artemis.app/Contents
java -mx2g -cp "$A/Java/artemis.jar:$A" uk.ac.sanger.artemis.components.ArtemisMain \
     S1.fasta + myproject_S1.busco.gff3.gz
```

Answer **Use index** to both questions, as above.

## What you see

- Each BUSCO hit is a `gene`, an `mRNA` and a `CDS`. The label is the BUSCO id.
- A CDS made of several miniprot segments is one feature with linked segments. The feature list
  at the bottom shows each segment on its own line. That is normal for GFF in Artemis.
- Click a CDS to see its qualifiers: `locus_tag`, `busco_id`, `busco_status`, `product`,
  `codon_start`, `translation_attached` and the same notes as in the GenBank file.
- The GFF3 has no `translation`. Artemis translates the CDS itself. The GenBank file keeps the
  checked `/translation`.
- Contigs with no BUSCO hit are still in the list. They show sequence only.
- Only BUSCO marker genes are annotated. Most genes in a genome are not BUSCO markers, so
  unannotated ORFs are expected.

## What does not work

These were tested with Artemis 18.2.0 and fail:

| Try | Result |
|---|---|
| Open the multi-record `.gbk` | Only the first contig is loaded; no contig list |
| Indexed FASTA + **uncompressed** GFF3 | Artemis tries to place every contig's genes on the chosen contig and stops with *"out of range location"* |
| Indexed FASTA + BGZF GFF3 **without `.tbi`** | No usable window |
| "Concatenate Sequences" instead of "Use index" | No usable window within 5 minutes on a ~46 Mb test assembly (and it joins the contigs into one artificial sequence) |

## Still available

- The main multi-record GenBank file (`<prefix>_<sample>.gbk`) is unchanged, for downstream tools.
- `--per-contig` still writes one GenBank file per contig, if you prefer GenBank in Artemis.

## Validation

`buscopainter-validate-artemis` re-reads the source FASTA, the painter GenBank and
`features.tsv`, and checks that:

- the FASTA copy is byte-identical and the `.fai` matches one recomputed (and `samtools faidx`,
  if installed);
- the GFF3 is valid BGZF, sorted, has a `##sequence-region` for every contig, unique IDs and
  existing Parents;
- every gene model has the same contig, gene span, strand and CDS segments as in the GenBank file
  and `features.tsv`, with the same BUSCO id, status and `codon_start`;
- each CDS segment's phase follows from `codon_start` and the segment lengths;
- the gene counts per contig agree in all three sources;
- the `.tbi` exists and `tabix` returns every contig's rows.

On the four validation assemblies (19,091 gene models, 1,849 of them multi-segment), every check
passed. The Artemis GUI itself was also driven by a script (`tools/artemis_gui_check/`). It
opened each assembly exactly as a user would, answered **Use index**, chose every contig in the
drop-down list, and scrolled the feature view along each contig. Every contig (725 in total) was
offered with the right length. Every gene model was loaded with the same coordinates, strand and
joins as in the GenBank file, and was drawn in the feature view. The menu route (File -> Open,
File -> Read An Entry) gave the same result.
