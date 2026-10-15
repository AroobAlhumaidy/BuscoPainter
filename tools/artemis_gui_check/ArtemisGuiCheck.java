// Drive the real Artemis GUI (not just its file reader) and record what a user would see.
//
// Starts Artemis exactly as `art <args>` does (ArtemisMain.main), answers the start-up dialogs
// with the button named by -Dclick (default "Use index"), waits for the Artemis window, then
// selects every entry of the contig drop-down (the real SequenceComboBox, so its own listeners
// run) and, for each contig, records:
//   contigs.tsv   : index, drop-down label, sequence length Artemis shows, features per entry
//   features.tsv  : every feature Artemis returns for the whole contig (key, location, locus_tag)
//   display.tsv   : features the feature view actually draws while it is scrolled along the
//                   contig window by window (the user's view), keyed by locus_tag
//   *.png         : the whole Artemis window painted offscreen for contigs named in <render>
//   dialogs.txt   : every dialog seen / button clicked
// Without a drop-down (e.g. a GenBank file) the single loaded sequence is reported as index -1.
//
// Usage: java -cp artemis.jar:<Artemis.app/Contents>:. [-Dclick="Use index"] ArtemisGuiCheck \
//            <outdir> <render contigs, comma-separated or -> <window bp> -- <art arguments...>
import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.ArtemisMain;
import uk.ac.sanger.artemis.components.EntryEdit;
import uk.ac.sanger.artemis.components.FeatureDisplay;
import uk.ac.sanger.artemis.components.SequenceComboBox;
import uk.ac.sanger.artemis.io.Range;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintWriter;
import java.lang.reflect.Method;
import java.util.*;
import java.util.List;

public class ArtemisGuiCheck {
  static final java.util.concurrent.ConcurrentLinkedQueue<String> chooserQueue = new java.util.concurrent.ConcurrentLinkedQueue<>();
  static final List<String> dialogs = Collections.synchronizedList(new ArrayList<>());

  public static void main(String[] a) throws Exception {
    File out = new File(a[0]);
    out.mkdirs();
    Set<String> render = new HashSet<>(Arrays.asList(a[1].split(",")));
    int window = Integer.parseInt(a[2]);  // 0 = one screen width
    String[] artArgs = Arrays.copyOfRange(a, 4, a.length);
    String click = System.getProperty("click", "Use index");
    int expectedEntries = Integer.getInteger("entries", 1 + (int) Arrays.stream(artArgs).filter("+"::equals).count());

    Thread clicker = new Thread(() -> {
      while (true) {
        try {
          Thread.sleep(150);
          SwingUtilities.invokeAndWait(() -> {
            for (Window w : Window.getWindows()) {
              if (!(w instanceof JDialog) || !w.isShowing()) continue;
              JFileChooser fc = findChooser(w);
              if (fc != null) {
                String next = chooserQueue.poll();
                if (next != null) {
                  dialogs.add("FILE CHOOSER '" + ((JDialog) w).getTitle() + "' -> " + next);
                  fc.setSelectedFile(new File(next));
                  fc.approveSelection();
                }
                continue;
              }
              String msg = describe((JDialog) w);
              AbstractButton b = null;
              String hit = null;
              for (String lab : click.split("\\|")) { b = findButton(w, lab); if (b != null) { hit = lab; break; } }
              if (b != null) { dialogs.add("CLICKED '" + hit + "' | " + msg); b.doClick(); }
              else if (!dialogs.contains("SEEN | " + msg)) dialogs.add("SEEN | " + msg);
            }
          });
        } catch (Exception e) { /* keep polling */ }
      }
    });
    clicker.setDaemon(true);
    clicker.start();

    long t0 = System.currentTimeMillis();
    String menuOpen = System.getProperty("menuOpen"), menuRead = System.getProperty("menuRead");
    if (menuOpen != null) {
      // user route: Artemis start window > File > Open ... (FASTA), then File > Read An Entry ... (GFF3)
      expectedEntries = menuRead == null ? 1 : 2;
      ArtemisMain.main(new String[0]);
      Frame start = null;
      while (start == null && System.currentTimeMillis() - t0 < 60_000) {
        for (Frame f : Frame.getFrames()) if (f instanceof ArtemisMain && f.isShowing()) start = f;
        Thread.sleep(200);
      }
      chooserQueue.add(menuOpen);
      clickMenu((JFrame) start, "Open ...");
      if (menuRead != null) {
        EntryEdit first = waitForEntryEdit(t0, 1);
        if (first == null) { System.out.println("NO_ARTEMIS_WINDOW"); writeDialogs(out); System.exit(4); }
        idle(); Thread.sleep(1000);
        chooserQueue.add(menuRead);
        clickMenu(first, "Read An Entry ...");
      }
    } else {
      ArtemisMain.main(artArgs);
    }
    EntryEdit ee = null;
    while (System.currentTimeMillis() - t0 < 300_000) {
      for (Frame f : Frame.getFrames())
        if (f instanceof EntryEdit && f.isShowing() && ((EntryEdit) f).getEntryGroup().size() >= expectedEntries)
          ee = (EntryEdit) f;
      if (ee != null) break;
      Thread.sleep(300);
    }
    if (ee == null) { System.out.println("NO_ARTEMIS_WINDOW"); writeDialogs(out); System.exit(4); }
    final EntryEdit frame = ee;
    idle(); Thread.sleep(1500); idle();
    long tLoad = System.currentTimeMillis() - t0;

    EntryGroup eg = frame.getEntryGroup();
    Method gfd = EntryEdit.class.getDeclaredMethod("getFeatureDisplay");
    gfd.setAccessible(true);
    FeatureDisplay fd = (FeatureDisplay) gfd.invoke(frame);
    Method visible = FeatureDisplay.class.getDeclaredMethod("getVisibleFeatures");
    visible.setAccessible(true);
    Method needUpd = FeatureDisplay.class.getDeclaredMethod("needVisibleFeatureVectorUpdate");
    needUpd.setAccessible(true);

    int win = window > 0 ? window : Math.max(1000, fd.getMaxVisibleBases() - 100);
    System.out.println("loaded in " + tLoad + " ms; view shows " + fd.getMaxVisibleBases() + " bases; sweep step " + win);
    SequenceComboBox combo = findCombo(frame);
    int n = combo == null ? 1 : combo.getItemCount();
    try (PrintWriter cw = new PrintWriter(new File(out, "contigs.tsv"));
         PrintWriter fw = new PrintWriter(new File(out, "features.tsv"));
         PrintWriter dw = new PrintWriter(new File(out, "display.tsv"))) {
      cw.println("index\tlabel\tartemis_length\tentries");
      fw.println("index\tlabel\tentry\tkey\tlocation\tlocus_tag\tID");
      dw.println("index\tlabel\tkey\tlocation\tlocus_tag");
      for (int i = 0; i < n; i++) {
        final int idx = i;
        if (combo != null) SwingUtilities.invokeAndWait(() -> combo.setSelectedIndex(idx));
        idle();
        String label = combo == null ? "-" : String.valueOf(combo.getItemAt(i));
        int len = eg.getSequenceLength();
        StringBuilder per = new StringBuilder();
        for (int e = 0; e < eg.size(); e++) {
          Entry en = eg.elementAt(e);
          FeatureVector fv = en.getFeaturesInRange(new Range(1, len));
          per.append(e == 0 ? "" : ";").append(en.getName()).append('=').append(fv.size());
          for (int k = 0; k < fv.size(); k++) {
            Feature f = fv.elementAt(k);
            fw.println(String.join("\t", String.valueOf(combo == null ? -1 : i), label, String.valueOf(en.getName()),
                f.getKey().toString(), f.getLocation().toString(), q(f, "locus_tag"), q(f, "ID")));
          }
        }
        cw.println((combo == null ? -1 : i) + "\t" + label + "\t" + len + "\t" + per);
        // sweep the feature view along the contig as a user scrolling would
        Set<String> seen = new LinkedHashSet<>();
        for (int centre = 1; centre <= len + win / 2; centre += win) {
          final int c = Math.min(centre, len);
          SwingUtilities.invokeAndWait(() -> fd.makeBaseVisible(c));
          SwingUtilities.invokeAndWait(() -> {
            try {
              needUpd.invoke(fd);
              FeatureVector v = (FeatureVector) visible.invoke(fd);
              for (int k = 0; k < v.size(); k++) {
                Feature f = v.elementAt(k);
                seen.add(f.getKey() + "\t" + f.getLocation() + "\t" + q(f, "locus_tag"));
              }
            } catch (Exception ex) { throw new RuntimeException(ex); }
          });
        }
        for (String s : seen) dw.println((combo == null ? -1 : i) + "\t" + label + "\t" + s);
        if (render.contains(label) || (combo == null && render.contains("-"))) paint(frame, fd, eg, new File(out, "window_" + label.replaceAll("[^A-Za-z0-9._-]", "_") + ".png"));
        if (i % 25 == 0) System.out.println("contig " + i + "/" + n + " " + label + " len=" + len + " " + per);
      }
    }
    writeDialogs(out);
    System.out.println("DONE contigs=" + n + " combo=" + (combo != null) + " load_ms=" + tLoad
        + " total_ms=" + (System.currentTimeMillis() - t0) + " visible_bases=" + fd.getMaxVisibleBases());
    System.exit(0);
  }

  static void paint(EntryEdit frame, FeatureDisplay fd, EntryGroup eg, File png) throws Exception {
    // centre on the first gene so the image shows annotation, then paint the real window
    FeatureVector all = eg.getAllFeatures();
    int centre = 1;
    for (int k = 0; k < all.size(); k++)
      if (all.elementAt(k).getKey().toString().equals("CDS")) { centre = all.elementAt(k).getFirstBase(); break; }
    final int c = centre;
    SwingUtilities.invokeAndWait(() -> fd.makeBaseVisible(c));
    idle();
    SwingUtilities.invokeAndWait(() -> {
      try {
        BufferedImage img = new BufferedImage(frame.getWidth(), frame.getHeight(), BufferedImage.TYPE_INT_RGB);
        Graphics2D g = img.createGraphics();
        frame.getRootPane().paint(g);
        g.dispose();
        ImageIO.write(img, "png", png);
      } catch (Exception ex) { throw new RuntimeException(ex); }
    });
  }

  static EntryEdit waitForEntryEdit(long t0, int entries) throws Exception {
    while (System.currentTimeMillis() - t0 < 300_000) {
      for (Frame f : Frame.getFrames())
        if (f instanceof EntryEdit && f.isShowing() && ((EntryEdit) f).getEntryGroup().size() >= entries) return (EntryEdit) f;
      Thread.sleep(300);
    }
    return null;
  }

  static void clickMenu(JFrame frame, String label) {
    JMenuItem it = null;
    JMenuBar bar = frame.getJMenuBar();
    for (int i = 0; bar != null && i < bar.getMenuCount() && it == null; i++) it = findItem(bar.getMenu(i), label);
    if (it == null) throw new RuntimeException("menu item not found: " + label);
    dialogs.add("MENU " + label);
    final JMenuItem item = it;
    SwingUtilities.invokeLater(item::doClick);   // later: the click opens a modal file chooser
  }

  static JMenuItem findItem(JMenu m, String label) {
    if (m == null) return null;
    for (Component c : m.getMenuComponents()) {
      if (c instanceof JMenu) { JMenuItem r = findItem((JMenu) c, label); if (r != null) return r; }
      else if (c instanceof JMenuItem && label.equals(((JMenuItem) c).getText())) return (JMenuItem) c;
    }
    return null;
  }

  static JFileChooser findChooser(Container c) {
    for (Component k : c.getComponents()) {
      if (k instanceof JFileChooser) return (JFileChooser) k;
      if (k instanceof Container) { JFileChooser r = findChooser((Container) k); if (r != null) return r; }
    }
    return null;
  }

  static String q(Feature f, String name) {
    try { String v = f.getValueOfQualifier(name); return v == null ? "" : v; } catch (Exception e) { return ""; }
  }

  static void idle() throws Exception {
    for (int i = 0; i < 3; i++) { SwingUtilities.invokeAndWait(() -> {}); Thread.sleep(40); }
  }

  static SequenceComboBox findCombo(Container c) {
    for (Component k : c.getComponents()) {
      if (k instanceof SequenceComboBox) return (SequenceComboBox) k;
      if (k instanceof Container) { SequenceComboBox r = findCombo((Container) k); if (r != null) return r; }
    }
    return null;
  }

  static AbstractButton findButton(Container c, String text) {
    for (Component k : c.getComponents()) {
      if (k instanceof AbstractButton && text.equals(((AbstractButton) k).getText())) return (AbstractButton) k;
      if (k instanceof Container) { AbstractButton r = findButton((Container) k, text); if (r != null) return r; }
    }
    return null;
  }

  static String describe(JDialog d) {
    StringBuilder sb = new StringBuilder(String.valueOf(d.getTitle()));
    collect(d, sb);
    return sb.toString().replace('\n', ' ');
  }

  static void collect(Container c, StringBuilder sb) {
    for (Component k : c.getComponents()) {
      if (k instanceof JLabel) sb.append(" | ").append(((JLabel) k).getText());
      if (k instanceof AbstractButton) sb.append(" [").append(((AbstractButton) k).getText()).append("]");
      if (k instanceof Container) collect((Container) k, sb);
    }
  }

  static void writeDialogs(File out) throws Exception {
    try (PrintWriter w = new PrintWriter(new File(out, "dialogs.txt"))) { synchronized (dialogs) { for (String s : dialogs) w.println(s); } }
  }
}
