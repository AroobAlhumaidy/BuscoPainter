# Artemis GUI check

Tests the whole-assembly Artemis view in the **real Artemis window**, not just its file reader.
It is not part of the package and needs a local Artemis install and Java 11+.

`ArtemisGuiCheck.java` starts Artemis as `art <files>` does, clicks **Use index** in the
start-up dialogs, chooses every contig in the drop-down list, scrolls the feature view along each
contig, and writes what Artemis loaded and drew (`contigs.tsv`, `features.tsv`, `display.tsv`,
`dialogs.txt`) plus screenshots of the window painted offscreen (`window_<contig>.png`).
A window opens on screen while it runs.

`compare_artemis_vs_genbank.py` compares that output with the painter GenBank file.

```bash
A=/Applications/Artemis.app/Contents          # or the folder that holds artemis.jar
javac -proc:none -cp "$A/Java/artemis.jar" ArtemisGuiCheck.java
E=results/myproject_S1_artemis
java -mx3g -cp "$A/Java/artemis.jar:$A:$PWD" ArtemisGuiCheck check_S1 contig_1 0 -- \
     "$PWD/$E/S1.fasta" + "$PWD/$E/myproject_S1.busco.gff3.gz"
python compare_artemis_vs_genbank.py results/myproject_S1.gbk check_S1     # PASS / FAIL
```

Arguments: `<output dir> <contigs to screenshot, comma-separated, or -> <scroll step in bp, 0 = one
screen> -- <files as for art>`. Use absolute paths.

To test the menu route instead (File -> Open ..., then File -> Read An Entry ...), leave the files
out after `--` and add `-DmenuOpen=<fasta> -DmenuRead=<gff3.gz>`.
`-Dclick="A|B"` sets which dialog buttons are clicked (default `Use index`).
