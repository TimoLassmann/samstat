# Code research: plots and reports (src/plot, src/report)

Scope: what SAMStat 2.3.0 shows today, how binning works, how output scales with read length, what a researcher sees first, and a keep/merge/drop/change verdict per plot.
Method: read src/plot/*, src/report/*, src/collect/collect.c (plot configuration), docs/json-format.md, docs/multiqc.md. Built Release in /private/tmp/claude-502/samstat-v24-plan/build-report. The dev/ files are 1-3 record fixtures, so I generated synthetic data with a script (/private/tmp/claude-502/samstat-v24-plan/gen.py): `short` (20k reads x 150 bp, ~1% error, Q~34) and `long` (1500 lognormal reads 0.5-30 kb plus one 200 kb read, 8% error, Q~14), each as SAM (with CIGAR + MD) and FASTQ, plus `mb.fastq` (one 1 Mb read + 2 short). I could not screenshot the HTML (browser tool blocks file://), so the visual findings come from the generated HTML/JS, not from rendering. Anything about how Plotly draws it is therefore inferred.

## 1. What the HTML report contains today

Order of sections (stat_report.c). Pairs repeat every section for Read 2.

| # | Section | Plot id | Plotly form | Data | Default visibility |
|---|---|---|---|---|---|
| 0 | Title, file size, **file mtime** (local time), "report produced by" | - | text | - | - |
| 1 | Mapping statistics (aligned only) / "File statistics" (FASTx: one line of read count) | MappingStatsTable | Plotly `table` (4 MAPQ bins + Total; count + %) | basic_nums | always |
| 2 | Length distribution | len_R1/R2 | line, **KDE** of counts (`tld_kde_count_pdf`), one series per MAPQ bin | per-length counts | all bins |
| 3 | Base composition | bc_R1/R2 | **stacked bars of raw counts** A/C/G/T(/N), one group of 4-5 series per MAPQ bin, a dropdown switches the bin | per position | first bin only (button to switch) |
| 4 | Base quality distribution | bq_R1/R2 | line: **mean Phred +/- standard error of the mean**, one line per MAPQ bin | per position | all bins |
| 5 | Error composition | mis_R1, ins_R1, del_R1 (+R2) | mismatch: stacked bars of raw counts by read base with MAPQ dropdown; ins and del: one line per MAPQ bin of raw counts | per position | mismatch first bin; ins/del all bins |
| 6 | Base comp / quality at read end (`--plotend`) | *_end | same as 3, 4 with x = distance from end, fixed 250 positions | | |
| 7 | Base comp / quality in 1% length intervals (`--norm`) | *_len_norm | same, 100 bins | | |
| 8 | Footer: About / Contact / Citation / Links | - | text | | |

Counts: 1 table + up to 5 plot types per read, 6 plots for single-end aligned (len, bc, bq, mis, ins, del), 12 for paired; 4 more per read with `--plotend`, 4 more with `--norm`. Unpaired FASTQ gives 3 plots (len, bc, bq); the table is replaced by a text line.

Facts that shape what the researcher sees:

* **First screen = file size and creation time, version string, then a numeric table of four MAPQ bins.** There is no headline (reads, mean length, N50, mean Q, % mapped, error rate), no pass/warn/fail, no per-plot interpretation hint beyond a one-line caption. The most distinctive SAMStat output (error profile by MAPQ) is the last block on the page.
* Bin labels are inconsistent: ids `mapq_ge30` -> label "MAPQ gt 30", `mapq_10_29` -> "MAPQ lt 30" (docs/json-format.md lists these verbatim); the mqc output uses ">= 30", "10-29".
* MAPQ bin 4 ("Unmapped") is silently dropped from every plot if empty (`plot_empty_series` is forced to 0 in `plot_data_config`, line ~968, so the caller's `PLOT_EMPTY_SERIES` flag is overwritten). For FASTQ the only series is named "Unmapped", which is a misleading label for plain sequencing reads.
* Mismatch/insertion/deletion plots show **raw counts**, not a rate per aligned base. Position-wise error rate (the quantity you compare between runs or against expectations) cannot be read off; coverage decay at read ends and read-length variation look like error trends. Only the base composition plot could be normalised but is also raw counts.
* Quality plot error bars are SEM (sd / sqrt(n)); with millions of bases they are invisible, so the plot is effectively a mean line. No median, quartiles or distribution (FastQC's box plots), no per-read mean-Q histogram. Mean of Phred values, not of error probabilities.
* HTML is not offline-capable: `<script src="https://cdn.plot.ly/plotly-1.58.5.min.js">` (html_header). Fixed `body{width:940px}`, `footer{position:absolute}` so mobile/narrow windows and long pages mis-render (inferred from CSS, not rendered).
* Report is not byte-deterministic (mtime + local time string in section 0); JSON and MultiQC outputs are deterministic.
* Text typos in captions ("Distibution", "compostion").

## 2. Binning (plot_interval.c)

Adapted from FastQC's BaseGroup. Applied to every position-indexed plot except the length KDE (`plot_add` skips `clu_data` for `PLOT_MOD_DENSITY`). `pd->len` = max read length + 1 over the whole file (`plot_data_resize_len`, collect.c lines 92-115, called for **all** plots of **both** read numbers).

| max len | rule | x axis |
|---|---|---|
| <= 250 | one point per position (labels 1..n) | numeric |
| 251-1000 | `get_plot_interval_lin`: first 250 single, then interval 5 x 10^k chosen so total bins < `target_n_clu`=500 (5, 10, 50, 100...) | categorical labels like "250-260" |
| > 1000 | `get_plot_interval_exp`: first 250 single, then step 10 (<2000), 50 (<5000), 100 (<10000), 500 (>=10000) | categorical |

Measured: 150 bp -> 155 points; 600-ish bp tail -> ~600 points; 200 kb read -> ~740 points per trace. So the *number of plotted points* is bounded (<~1000) and the HTML stays small (long.sam: 320 KB, 1500 reads incl. one 200 kb; mb.fastq with a 1 Mb read: 109 KB).

Problems:

1. **Binned values are `round(mean count per position)`** (`summarise_data`): sums of counts divided by number of positions in the bin, then rounded to an integer. Sparse tails (positions covered by few reads; and each of the 42 quality rows separately) round to 0 or 1, distorting the quality mean and dropping tail composition. By inspection of the code; not quantified. Averaging also means bar heights in stacked bars are per-position averages, not totals, so they are comparable only inside one plot.
2. The long-read axis is a *position in read* axis that mixes extreme coverage change (a 5 kb median read with a 200 kb outlier: after ~5 kb only the outlier contributes). The initial `range: [0,250]` (plot_add_code; set when any series has > 250 points) means a long-read report **opens zoomed on the first 250 bases**, with a rangeslider to scroll. For long reads that first window is the least informative part of the data (adapter/start effects only).
3. Binning is only for position-axis plots, not for length: the KDE is computed on every distinct length (`tld_kde_count_pdf`, x labels printed with `%f`), unbounded in a pathological case (all-distinct lengths: points = number of distinct lengths, up to 1.5k here; 100k+ distinct lengths in a big run would be heavy). The KDE smooths away sharp length peaks (e.g. amplicons, 16S, adapter-dimers) and has no log scale; the x-axis is linear so a 100 bp - 100 kb distribution is unreadable.
4. Three nearly identical binning code paths (lin/exp/none) plus a hard-coded 250 in collect.c (`plot_data_alloc(...,250,...)`) and in plot.c; a lot of commented-out dead code (~600 lines in plot.c).

## 3. How output and memory scale with read length

Measured (Release, this machine):

| input | wall | peak RSS | HTML | JSON |
|---|---|---|---|---|
| short.sam 20k x 150 | 0.05 s | 28 MB | 63 KB | 78 KB |
| short.fastq | 0.02 s | 18 MB | 15 KB | 27 KB |
| long.sam (1500 reads + one 200 kb) | 1.9 s | 877 MB | 321 KB | 46 MB |
| long.fastq | 0.7 s | 814 MB | 99 KB | 20 MB |
| mb.fastq (reads of 1 Mb, 5 kb, 3 kb) | 3.6 s | **3.8 GB** | 109 KB | **96 MB** |
| long.sam `-p 100` (peek) | 0.9 s | 835 MB | | |

* **Memory is O(max read length), independent of the number of reads.** Dense `uint64_t[L][len]` matrices: per read number about 5x4 (bc) + 42x4 (qual) + 5x4 (mis) + 4 + 4 + 4 rows = 220 rows x 8 B = ~1.8 KB per position; R1 and R2 both allocated -> ~3.5 KB/position -> 3.5 GB for a 1 Mb read (matches 3.8 GB measured). One ultra-long read (ONT up to several Mb) can exhaust a laptop. `-p` does not help because the matrices grow to the longest read seen.
* **HTML size is bounded** (binning), good. **JSON size is not**: 96 MB for one 1 Mb read (quality matrix 42 rows per position, trimmed only at the tail), and this is the format the MultiQC/automation story is built on. MultiQC `_mqc` files are bounded (2000-point blocks, 25-35 KB).
* Speed: ~100k short reads/s-class on one core for the parse path (20k reads 0.05 s); for the long-read case the time is dominated by CIGAR/MD walk and matrix allocation/finalisation (3.6 s for 1 Mb read, almost all allocation and JSON output). Not measured on a real multi-GB BAM; UNVERIFIED how the htslib read path compares to samtools.
* JSON `length` is one array of `max_len` integers per MAPQ bin (so 1 Mb read -> 4 arrays of up to 1M entries, mostly zeros).

## 4. JSON and MultiQC outputs (src/report/json_report.c, mqc_report.c)

* JSON (format_version 1): the raw dense counts (length hist, base comp, quality 42xposition, mismatch/ins/del per position) per read number and MAPQ bin. Deterministic, exact. Good for tests and downstream tools, but it is a *dump of the internal matrices*, not a report: no derived metrics (mean length, N50, mean Q, error rate), and size is O(max len) as above. It is the only place the full-resolution data lives (HTML is binned, rounded).
* MultiQC custom content: 4 files per input (general stats columns: reads, % MAPQ>=30, % unmapped, % properly paired; MAPQ table; read-length linegraph; mean-Q-per-position linegraph, read 1 only, all MAPQ bins merged). **The error profiles (SAMStat's distinguishing feature) are not exported to MultiQC at all**, nor base composition, nor read 2, nor per-MAPQ splits. In a multi-sample run via MultiQC SAMStat adds a length and a quality line only, similar to what FastQC/fastp modules already show. The length series is summed in 2000-point blocks (x = first value of block), which for a 1 Mb read means 500 bp blocks; fine.
* Both outputs are opt-in (`--json`, `--multiqc`), defaults off. No native MultiQC module (needs a parseable JSON/TSV that MultiQC's module parser accepts; UNVERIFIED whether the existing custom-content route is sufficient for the nf-core module; see https://docs.seqera.io/multiqc/custom_content for the format used).

## 5. Verdict per plot / output

Guiding principle: a fixed, small, FastQC-like set of plots, each readable at a glance. Proposed target set for v2.4 (mapped or unmapped, short or long): header summary + 6 plots (length, quality, base content, and the three error profiles merged into one, plus mapping/MAPQ). Verdicts:

| Current item | Verdict | Reason / change |
|---|---|---|
| Header (file size, mtime, version) | **change** | Replace with a one-screen summary: reads, bases, mean/median length, N50 (new), mean Q (probability-averaged), % mapped, % MAPQ>=30, read identity/error rate if MD or `NM` present; flag pass/warn/fail like FastQC? (maintainer decision; keep it to a few badges at most). Drop mtime to restore determinism. |
| MAPQ table (Plotly table) | **change** | Keep the content (4 bins x R1/R2) but render as a plain HTML table or a single horizontal stacked bar; no Plotly dependency; fix bin labels (">= 30", "10-29", "< 10", "Unmapped"); paired info as text. For FASTQ, drop the "Unmapped" label. |
| Length distribution (KDE) | **change** | Replace KDE with a real histogram, log-x for long reads (bin automatically: linear <= ~1 kb range, log otherwise), mark N50. Per-MAPQ split only if aligned. Biggest value-add for long reads. |
| Base composition (stacked counts, dropdown by MAPQ) | **change** | Show A/C/G/T/N as **% per position**, one view (all reads), drop dropdown and per-MAPQ groups (or keep MAPQ split only for error plots). Bin by position intervals as now but with exact (non-rounded) aggregation. Likely little value for long reads beyond start/end: fine at coarse bins. |
| Base quality (mean +/- SEM per MAPQ) | **change** | Show mean and a spread (median + IQR or 10-90%) per binned position, like FastQC's box plot; add **read-level mean-Q histogram** (long reads: the standard QC, e.g. NanoPlot/pycoQC) as the second panel. Drop SEM. |
| Mismatch plot (stacked counts by read base, MAPQ dropdown) | **merge + change** | Merge mismatch + insertion + deletion into **one "Error profile along the read" plot**: three lines (mismatch, ins, del) as **rate per aligned base per position**, with a MAPQ >=30 vs <30 split as the only stratification. This is the distinctive plot; it should be near the top, not last. Per-base substitution type (A->C etc.) can be a small 4x4 heatmap/table if the maintainer wants to keep substitution detail, else drop. |
| Insertion plot | **merge** | Into the error profile above. |
| Deletion plot | **merge** | Into the error profile above. |
| Read 2 duplicate sections | **change** | Overlay R1/R2 as two lines in the same plot (or a toggle) instead of doubling the page (12 plots -> 6). |
| `--plotend` plots (4) | **drop** | Rarely needed; for long reads the end-aligned view is the relevant one, so replace by showing both start and end windows in the position plots (e.g. first 100 and last 100 bases) rather than a separate flag. |
| `--norm` plots (4) | **drop** or **make the long-read default** | 1%-of-length axis is the right axis for variable-length long reads (position in read is meaningless past a few kb). If kept, make it the automatic axis when the length distribution is wide; delete the flag. |
| Partial-report red warning | **keep** | Make it a banner at the top, not in the middle. |
| Footer | **change** | Shrink to a citation line; remove `position:absolute` footer. |
| Plotly CDN | **change** | Embed minimal plotting (hand-rolled SVG, a few KB) or inline a trimmed Plotly basic bundle (~1 MB) so the report works offline and is self-contained; HTML size currently 15-320 KB so a 1 MB inline bundle would dominate. Maintainer call; see section 6. |
| JSON (`--json`) | **change** | Make it small and derived (summary metrics + the same *binned* arrays as the plots, bounded size) so it is the machine twin of the HTML; the exact full-resolution dump is not needed by researchers. This also gives MultiQC a good source. Breaks `format_version 1` -> bump to 2. |
| MultiQC custom content (`--multiqc`) | **change** | Add the error-profile linegraph, per-read-mean-Q or read-length with N50; consider writing the `_mqc` files by default into the same file set (cheap) or merge into JSON. Do not add more sections than the HTML has. |

## 6. Key recommendations for the planners

1. **Fix the memory model before anything else for long reads**: dense `[L][maxlen]` matrices x R1/R2 -> 3.8 GB at 1 Mb. Aggregate into bounded bins while parsing (e.g. a fixed ~400 position bins: exact first 100-250 positions + log-spaced after, or 1%-of-length bins) and keep a per-read histogram for length and mean Q. Output size then stops depending on the longest read, and the rounding bug in `summarise_data` disappears because aggregation happens on exact counters.
2. **Rate, not counts, for error plots**, per aligned base, with a MAPQ split. This is the one thing no other QC tool gives (see the paper) and is currently the least readable part of the page.
3. **Make the first screen the answer**: a single summary block (reads, bases, N50/mean length, mean Q, % mapped, % MAPQ>=30, error rate) + length histogram + error profile. Everything else below.
4. **One page layout for short and long reads**; the tool picks axis scale from the data (position bins for short reads; length-normalised or log bins for long), no flags.
5. **Offline, deterministic, one HTML file**; keep Plotly only if a trimmed build is acceptable, otherwise tiny hand-written SVG (6 plots of lines/bars/histogram do not need Plotly; also kills the 1.58.5 CDN pin).
6. **Delete ~600 lines of commented dead code** in plot.c, the three duplicated binning paths, and the R1/R2 duplication when touching these files.
7. Keep stable outputs: JSON deterministic, exit codes, file names, `--json`, `--multiqc` (CI and bioconda tests depend on them; goldens in tests/golden and tests/json/golden will all change with these plot changes).

## 7. Unverified / not done

* No rendered screenshots (browser tool blocked file:// URLs); layout claims come from CSS/JS source.
* No test on a real large BAM or real ONT/PacBio data; numbers are from synthetic reads with random sequence.
* Wall times are single runs on one machine.
* No external tool feature claims beyond the MultiQC custom-content doc link cited above (not re-verified online).
