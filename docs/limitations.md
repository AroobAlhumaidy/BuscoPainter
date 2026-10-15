# Known limitations

## Scientific

1. **BuscoPainter transfers BUSCO annotations; it is not a gene predictor.** Features mark where
   BUSCO's gene predictor (miniprot) aligned BUSCO marker proteins. Only a few thousand conserved
   genes are covered. The structures, ends and products are not curated. Most open reading frames
   in a genome get no feature, and that is expected: in each validation assembly, thousands of
   long ORFs do not overlap any BUSCO feature.

1a. **BUSCO reports only its scored markers, and the choice of lineage matters.** A lineage's
   `refseq_db.faa` contains more ortholog groups than it has HMMs (euglenozoa_odb12: 3,991 vs
   3,611). miniprot aligns all of them, but BUSCO scores and reports only the 3,611 with HMMs, so
   alignments to the other 380 groups never reach `full_table.tsv` or BuscoPainter. A gene can
   therefore be a marker in one lineage (for example `leishmaniinae_odb12`) and absent from
   another. Repeat-rich proteins (WD40, leucine-rich repeat) are often reported as Fragmented or
   Missing, or as a gene model that only partly covers the open reading frame.

2. **Multi-segment models are unverified.** About 10 % of models in the validation data join ≥ 2 miniprot CDS segments, often across gaps > 1 kb. In organisms with few
   cis-introns, such as trypanosomatids, these are more likely to be frameshifts, assembly errors or
   spurious spliced alignments than real introns. They are exported exactly as in the GFF and
   flagged `MULTI_SEGMENT_MODEL`. Treat them with caution in any gene-structure analysis.

2a. **BUSCO-trimmed models.** BUSCO 6 sometimes drops terminal segments of a miniprot model: its
   `full_table.tsv` span and `.faa` protein then cover only the remaining whole CDS rows, while the
   per-BUSCO GFF still holds the whole model. Since T7.1-fixed.3 only BUSCO's rows are exported
   (`FULLTABLE_SPAN_TRIMMED`; under 1 % of models in the validation data). If the span does not fall on
   whole rows, the whole GFF model is exported without a translation
   (`FULLTABLE_SPAN_DIFFERS_FROM_GFF`; none in the validation data). When the first kept row has a
   GFF phase > 0, BUSCO's `.faa` protein has one extra N-terminal residue (the codon that spans the
   dropped junction); the GenBank translation starts at `/codon_start` and omits it.

3. **Frame-inconsistent models** (the GFF phases do not form one reading frame) and models with
   miniprot `Frameshift` or `StopCodon` attributes are exported without `/translation`.

4. **Gene ends are not asserted.** Many models do not start with ATG, and
   some have no stop codon, because miniprot aligns protein homology, not complete ORFs. Locations
   carry no `<` or `>` partial markers. `starts_with_ATG` and `ends_with_stop` in `features.tsv`
   record what was observed (see output_formats.md for the caveat on `ends_with_stop`).

5. **BUSCO's internal-stop-codon statistic ("E") is not reproduced.** In the validation data, BUSCO 6
   reported E for several hundred BUSCOs per sample. This number could not be derived from BUSCO's
   selected proteins, miniprot `##STA`, or PAF tags. BuscoPainter neither uses nor changes it. Its
   own check found in-frame stops in far fewer exported models, and these were exported without
   translation.

6. **Standard genetic code only.** `/transl_table=1` is fixed. For organisms with an alternative
   nuclear code (for example, some ciliates and yeasts), translations would be wrong or omitted.
   Do not use those translations.

7. **Product names** are the BUSCO `full_table.tsv` `Description`, copied verbatim (for example,
   `hypothetical_protein_-_conserved`). They are not functional annotations.

8. **Fixed taxon wording.** The multi-segment note always says "cis-introns are rare in
   trypanosomatids", whatever the organism.

## Compatibility

9. **Tested only with BUSCO 6.0.0 + miniprot 0.18-r281, genome mode, `euglenozoa_odb12`.**
   MetaEuk or AUGUSTUS GFFs, BUSCO 5.x, other lineages, and protein or transcriptome modes have not
   been tested.

10. **The FASTA must be the BUSCO input.** Matching uses contig names only. A different assembly
    with the same names would be painted wrongly without any error. Compare checksums.

11. **No whitespace in paths.** Config, FASTA and BUSCO paths with spaces or tabs are rejected.
    Relative paths are resolved from the current working directory.

12. **Long contig IDs.** IDs longer than 16 characters give a non-standard GenBank `LOCUS` line. IDs
    containing `.` are split into ACCESSION/VERSION by Biopython. The file still round-trips in
    Biopython and passes the validator, but other readers (including Artemis) have not been tested
    with such IDs.

13. **Exit code 0 is not proof of success.** A failed count reconciliation gives exit code 3, but
    rejected hits are reported, not fatal. Always read the log and `skipped.tsv`, and run the
    validator.

14. **OrthoDB API lookup** (`use_orthodb_api`) was not used in the validation. Per-contig export
    was validated in T7.1-fixed.2 and T7.1-fixed.3 (every contig file read with Biopython and with
    Artemis 18.2.0's file reader; the Artemis GUI itself was not driven by a test). The
    whole-assembly export (`buscopainter-artemis`) was tested in the Artemis 18.2.0 GUI on macOS
    only, driven by a script rather than by hand. It needs `tabix`; without the `.tbi` index Artemis
    cannot use the GFF3 (docs/artemis.md). Other Artemis versions and IGV/JBrowse were not tested.

15. **The validator** reads JSON configs only, needs an uncompressed `full_table.tsv`, and needs the
    `busco_sequences/` layout.

16. **Platform.** Validated on macOS (arm64) only. The code is pure Python and should run elsewhere,
    but that has not been tested.
