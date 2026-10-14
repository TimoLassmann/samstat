# SAMStat 2.4: presentation (final design)

Base: the "researcher" proposal. Grafted from the "minimal" proposal: the fixed layout that never shifts (plots that do not apply become a one-line note), a single small set of SVG primitives, and a smaller check table. Cut from both: an `index.html` of links, the fraction-by-bin quality variant, the stacked MAPQ bar, a footer "thresholds link", and the per-plot data tables beyond one `<details>` block.

## 0. What the report is for

A researcher reads the report for 30 seconds and must be able to answer four questions:

1. Did the run work?
2. Is this short-read or long-read data?
3. How well did the reads map?
4. Where along the read do errors occur, and do they differ between well and poorly mapped reads?

Every element below serves one of these. There is one page and one order, and nothing about the layout can be configured.

## 1. Code facts this design relies on (verified by reading the source)

- `src/plot/plot.c:1414` emits `<script src="https://cdn.plot.ly/plotly-1.58.5.min.js">`. This is the only network dependency the report has.
- `src/plot/plot.c:243-248` builds a Plotly `updatemenus` dropdown, which is the base-composition dropdown.
- Report size today: `plot.c` has 1,635 lines, `plot_group.c` 263, `plot_interval.c` 333, `plot_test.c` 105 and `stat_report.c` 659.
- The MAPQ table is a Plotly table, rendered by `<script>` at `stat_report.c:487`.
- Output name is `<input>.samstat.html`, or `<outdir>/<basename>.samstat.html` (`stat_report.c:263-265`). This does not change.
- Base quality is capped at Q41 (`collect.c:43`, `QUAL_MAX_BIN 41`).
- The only verbosity option is `--verbose` (`param.c:53`). There is no `--quiet`, and none is added.

## 2. Decisions

### D1. One self-contained HTML file with inline static SVG generated in C

There is no Plotly, no CDN and no JavaScript.

- **Why:** a CDN report fails on cluster nodes, in containers, offline, and in archives.
  - Inlining Plotly instead costs 1.2 MB (basic bundle), 1.51 MB (cartesian) or 4.83 MB (full), per https://github.com/plotly/plotly.js/blob/master/dist/README.md (verified in `research/presentation.md`).
  - NanoPlot's offline mode adds about 4.8 MB per file (`research/nanoplot-cramino.md`).
- **What SAMStat needs:** only line plots, step histograms and simple bars, each with at most a few hundred points after binning.
  - Each plot is about 5-20 KB.
  - Hover values come from SVG `<title>` elements.
  - The page prints as vector graphics.
- **Primitives:** one new file, `src/plot/svg_plot.{c,h}`, of about 400 lines, containing:
  - a nice-ticks helper
  - linear and log axes
  - a polyline with dash pattern
  - step bars for histograms
  - direct end-of-line labels
  - `<title>` tooltips
- **Not built:** scatter, heatmap, box, violin, zoom, pan and trace toggles.
- **Fallback:** inline plotly-basic, which keeps the rest of this design unchanged.

### D2. Fixed page structure, identical for every input

1. **Header line:** file name, detected format, read type and SAMStat version.
   - Read type is "short" or "long", auto-detected (long when median length is about 1 kb or more; the rule is owned by the collector topic).
2. **Summary strip (numbers only):**
   - reads, bases, median length, N50
   - mean Q, computed by averaging error probabilities
   - mapped %, MAPQ≥30 %, median identity
   - For FASTA/FASTQ the aligned-only fields are omitted, not shown as 0.
3. **Checks table:** at most 5 rows (D4).
4. **Plots, in this order:**

   | # | Title (a question) | Content | Input |
   |---|---|---|---|
   | 1 | How long are the reads? | Log-x histogram with N50 marked, and two traces: reads and bases | all |
   | 2 | How good are the reads overall? | Per-read mean-Q histogram | all |
   | 3 | How well did reads map? | MAPQ bins as a plain HTML table (counts, %) | aligned |
   | 4 | Where along the read do errors occur? | Mismatch, insertion and deletion **rates per aligned base**: three stacked panels sharing an x-axis, one line per MAPQ bin. Given the most vertical space. | aligned |
   | 5 | How accurate are the aligned reads? | Per-read identity histogram, one line per MAPQ bin | aligned |
   | 6 | Does quality drop along the read? | Median and IQR of base quality by position | all |
   | 7 | Is base composition even along the read? | % A/C/G/T/N by position. No dropdown. One small-multiple panel per MAPQ bin. | all |

   - **Plot 4** is SAMStat's differentiator. It changes from counts to rates, which answers issue #12.
   - **Plot 7** stays because it is the paper's lead example (an N spike at position 40 in poorly mapped reads, `research/paper-samstat.md`). It is placed last. See DECISION NEEDED 2.
   - **R1/R2:** these are two traces within the same panel (R1 solid, R2 dashed), never separate plots. R2 is not drawn for single-end input.
   - **Plots that don't apply:**
     - Unaligned input: plots 3-5 are replaced by a one-line note ("not applicable: input is not aligned").
     - Aligned input without MD: plots 4-5 are replaced by "not available: no MD tag (re-run `samtools calmd`)". If the collector topic adds an MD-free path, this note disappears.
   - Either way the order and numbering never change.
5. **Footer:** SAMStat version, full command line, and the citation.

### D3. Help text for every plot

Each plot has three pieces of text:

- a title phrased as a question;
- a one-line caption giving the typical range, never a pass/fail statement;
- a collapsed `<details>` "How to read this", three lines long: what good looks like, what bad looks like, and the common benign cause.

Under that, one `<details>` "Data" block holds the binned numbers as an HTML table, for copy-paste and screen readers.

All titles, captions, help text and the check thresholds live in one file, `src/report/report_help.c`, so the wording can be reviewed in one place.

Typical-range numbers (for example Illumina vs ONT mismatch rates) are UNVERIFIED until they are measured on real NovaSeq, Dorado and Revio files. Until a range is measured, the caption describes shape only.

### D4. Few advisory checks

A check is shown as a shape plus a word (OK ✓, CHECK ▲, FAIL ■), followed by the measured value and one sentence. Colour is never the only cue.

| Check | Rule (initial, to calibrate) |
|---|---|
| Input read completely | FAIL on empty, truncated or unterminated input |
| Mapped fraction (aligned only) | CHECK below 50 %, FAIL below 10 % |
| MAPQ<10 fraction (aligned only) | CHECK above 30 % |
| Median read Q | CHECK below a platform-aware floor (short Q20, long Q10). UNVERIFIED. |
| Median identity (aligned, MD present) | CHECK below a floor. Value UNVERIFIED. |

- **Limit:** 5 checks for v2.4. Adding another needs a written justification.
- **Info only, never flagged:** GC, duplication, base composition, length spread, secondary/supplementary fraction, modification calls. These are FastQC's false-alarm sources (`research/fastqc.md`).
- **Where thresholds live:** one `static const` struct in `report_help.c`, mirrored in a table in `docs/reading-the-report.md`.
- **No FAIL where unsure:** a platform-dependent check never goes beyond CHECK until it has been calibrated.
- **Exit codes:** these stay 0/1/2 and keep their I/O and usage meaning. Checks never change the exit code.

### D5. Terminal summary line

One line per input is printed to stderr, always, with no flag:

```
sample.bam: 12.4M reads, N50 18.2 kb, mean Q 17.1, 94.2% mapped, identity 96.8% [OK] -> sample.bam.samstat.html
```

- **Field order:** the same as the summary strip. Aligned-only fields are omitted for FASTQ.
- **Status field:** `[OK]`, or `[CHECK: <first check name>]`, or `[FAIL: <first check name>]`.
- **Formatting:** plain text with no ANSI colour.
- **Progress:** a progress line appears only when stderr is a TTY (`research/code-cli.md`, item 6).

### D6. Palette and accessibility

- **Okabe-Ito palette:** https://jfly.uni-koeln.de/color/ (page not fetched, so the citation is UNVERIFIED; the palette is widely documented).
- **MAPQ colours and dashes:** each MAPQ bin gets one fixed colour and dash pattern, the same on every plot.
  - ≥30 solid, 10-29 dashed, <10 dotted, unmapped grey dash-dot.
  - This keeps grayscale printouts readable.
- **Labels:** lines are labelled directly at their end, and a legend is used only where lines overlap.
- **Text:** system font stack (no web fonts), at least 12 px, and contrast of at least 4.5:1 (https://www.w3.org/TR/WCAG21/#contrast-minimum).
- **SVG markup:** every SVG has `role="img"`, `<title>` and `<desc>`.
- **Theme and width:** light theme only, and the layout works at a 360 px width.
- **Print:** one `@media print` block that forces `<details>` open, sets `break-inside: avoid` on figures and targets A4 portrait. There is no PDF generator.

### D7. Multi-sample comparison is left to MultiQC

- **No comparison inside SAMStat:** no comparison view and no `index.html` of links.
- **MultiQC:** `--json` and `--multiqc` stay the multi-sample path. MultiQC's General Statistics table is the multi-sample summary (https://docs.seqera.io/multiqc/reports/customisation, verified in the researcher proposal).
- **Content to add (MultiQC topic):** the binned error profiles by MAPQ, plus the length and identity summaries. The 2.3 `_mqc` output omits the error profiles (`research/code-report.md`).
- **One model:** HTML, JSON and the terminal line all render from the same in-memory stats model, so their headline numbers cannot disagree.

## 3. Explicitly OUT

- Plotly or any other charting JS, any CDN, and web fonts.
- Interactive features: dropdowns, tabs, zoom/pan, trace toggles, a theme switcher and dark mode.
- `--plotend` and `--norm` as presentation switches. The axis is chosen automatically by the collector.
- Configurable thresholds, colours, bin sizes or MAPQ bins.
- Per-tile, k-mer, overrepresented-sequence, adapter and GC plots or flags.
  - Adapter content is a separate scope decision, not part of this topic.
- 2D length-vs-quality heatmaps, KDE and hex scatter.
- The fraction-by-bin quality variant for binned (NovaSeq X) data. It is deferred because it adds a second form of the same plot.
- An in-report comparison view or an `index.html`.
- A whole-file pass/fail verdict.
- A `--quiet` flag.

## 4. Dependencies

- **Collector topic (blocking):**
  - It must produce bounded position bins, with counts aggregated while parsing.
  - It must provide the log-length histogram with bases per bin, the per-read mean-Q histogram, and the identity histogram by MAPQ.
  - It must raise the Q41 cap.
  - Without bounded bins the SVG size is unbounded. The collector work and this work must land together.
- **JSON/MultiQC topic:** the new summary fields are shared. The JSON format change is owned there.
- **Order of work:** collector, then JSON, then HTML and the terminal line.

## 5. Affected files

| File | Change |
|---|---|
| `src/plot/svg_plot.{c,h}` | **New.** The SVG primitives (D1). |
| `src/report/report_help.c` | **New.** Titles, captions, help text and the check threshold table. |
| `src/plot/plot.c`, `plot_group.c`, `plot_test.c` | Replace the Plotly emitter with SVG calls. Most of `plot.c` is deleted. |
| `src/plot/plot_interval.c` | Keep only what still fits once bins come from the collector (`clu_data` late binning goes). |
| `src/report/stat_report.c` | New layout: header, summary strip, check table, fixed plot order, HTML MAPQ table, footer, print CSS. |
| `src/samstat.c` | Terminal summary line, and progress output when stderr is a TTY. |
| `docs/reading-the-report.md` | **New.** One page explaining each plot, plus the thresholds table. |
| `tests/golden`, `tests/report` | Regenerate the goldens and add the new tests. |
| `json_report.c`, `mqc_report.c` | No layout change. Only shared summary fields are added (owned by the JSON/MultiQC topic). |

## 6. Effort (one maintainer)

| Item | Days |
|---|---|
| SVG primitives and tick logic | 2-3 |
| Layout, CSS, print styles, help text | 2 |
| Wiring the 7 plots and the summary strip | 2 |
| Checks and terminal line | 1 |
| Tests and one golden regeneration | 1-2 |
| **Total** | **8-10**, on top of the collector work |

## 7. Risks

- **Hand-written SVG layout** (tick labels, label overlap):
  - Keep the primitive set tiny and add golden tests on the SVG text.
  - If the work overruns, switch to inline plotly-basic.
- **Rendering was never checked visually** during research.
  - Before release, do a manual check in Chrome and Firefox plus one print-to-PDF.
- **Wrong long-read thresholds would destroy trust**, as FastQC's did.
  - Every platform-dependent check stays at CHECK at most.
  - Calibrate on real ONT, HiFi and NovaSeq X data before release.
- **Interactivity is lost:** there is no zoom, and `<title>` tooltips are plainer than Plotly's.
  - This is accepted for QC use. That zoom is rarely used is UNVERIFIED.
- **Binned qualities:** median/IQR on 3-bin NovaSeq X data looks blocky. This is accepted for 2.4.
- **MultiQC consumers:** extend the `_mqc` fields, never rename them.

## 8. Tests

- **Fixed structure:** every report (short BAM, long BAM, FASTQ, BAM without MD) contains the 7 plot titles in order, with the n/a notes where they apply.
- **Offline:** no `<script src`, and no `http(s)://` inside `src=` or `href=` attributes for resources. Text links to docs and the citation are allowed.
- **SVG:**
  - Each SVG parses as XML (a CMake script check).
  - Each SVG has `role="img"`, `<title>` and `<desc>`.
  - Each series stays under its point cap.
- **Size bound:** a report for a fixture with a 1 Mb read stays under 500 KB. This replaces the 4 MB golden skip.
- **Checks:** a unit test feeds synthetic stats at each threshold boundary and asserts the word, shape and value.
- **Terminal line:** stderr goldens for the OK, CHECK and FAIL cases, and for a FASTQ input.
- **Consistency:** the HTML summary strip, the JSON and the terminal line give identical headline numbers for the same fixture.
- **Print:** the report contains an `@media print` block.
- **Goldens:** regenerated in one dedicated commit.

## 9. DECISION NEEDED

1. **SVG or Plotly:** accept the hand-written SVG emitter (recommended), or use inline plotly-basic (+1.2 MB per report)?
2. **Base composition:** keep it as plot 7 (recommended, since it is the paper's lead example), or cut it by default?
3. **Check thresholds:** sign off the initial values in D4, including the platform-aware median-Q floors, before calibration data is available. Or ship those two rows as info-only until they are calibrated.
4. **MAPQ split for plot 7:** show base composition split by MAPQ bin (small multiples, as in the paper), or as a single panel over all reads (simpler)?
