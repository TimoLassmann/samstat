# SAMStat 2.4: simplicity review (cross-topic)

Reviewer role: simplicity gatekeeper. Inputs: all nine files in `docs/plan-v2.4/topics/`. No external facts were added in this review; every URL and every UNVERIFIED tag stays as the topic files have it.

The rule applied throughout: **a FastQC user should be able to read the whole report in 30 seconds the first time they see it.** Every element must answer one of four questions: did the run work, short or long reads, how well did reads map, and where do errors occur by mapping quality. Anything else is cut or deferred.

Where this review contradicts a topic file, this review wins for the implementation workflow. The maintainer still signs off the items in section 7.

---

## 1. Surface budget (firm)

| Surface | 2.3.0 | 2.4 budget | Notes |
|---|---|---|---|
| Plots in the HTML | 9-10 Plotly figures, with a dropdown | **6, fixed order** | Plus one 4-row MAPQ table. Unaligned input gets plots 1-4. |
| Summary strip values | file size/mtime header | **8** | The aligned-only values are omitted for FASTA/FASTQ. |
| Status flags (OK/CHECK/FAIL) | 0 | **0** | Factual notes only (section 4). |
| Documented CLI options | 12 (`-d -p -t -l --seed --norm --plotend --json --multiqc --verbose -h -v`) | **8** (`-d -p -t --json --multiqc --verbose -h -v`) | Plus `-` for stdin. Four deprecated options are accepted as hidden no-ops. **Zero new options.** |
| Output files per input | html, optional json, optional `_mqc` set | **unchanged names** | `stdin.samstat.*` is the only new naming rule. |
| MultiQC general-stats columns | 4 | **7** | Adds `samstat_n50`, `samstat_mean_q` and `samstat_median_identity`. |
| MultiQC sections | 4 | **5** | Adds the `samstat_error` linegraph. |
| JSON format | v1 | **v2, one break, then frozen for 2.x** | |
| Stderr output per input | varies | **1 final line**, plus a progress counter only on a TTY | |

Adding anything beyond this budget needs a written justification in the PR, the same rule presentation set for checks.

---

## 2. The standard plot set (final)

Page order: header line, summary strip, notes (only if any), MAPQ table (aligned only), then plots 1-6, then the footer.

**Summary strip (8 values):** reads, bases, N50, median length, mean Q (error-probability averaged), % mapped, % MAPQ >= 30, median identity. Secondary and supplementary counts go in a footnote as raw counts.

**MAPQ table, "How well did reads map?":** it has 4 rows (>=30, 10-29, <10, unmapped). Each row shows the line style used in every plot (colour plus dash), the read count and the %. It is the legend for every MAPQ line on the page, which is why it comes before the plots. It replaces standard-plots' stacked bar (plot 5) and is presentation's HTML table.

| # | Title (a question) | Content | Input | FastQC equivalent |
|---|---|---|---|---|
| 1 | How long are the reads? | Log-x length histogram (reads), with an N50 line. One line per MAPQ class. | all | Sequence Length Distribution |
| 2 | How good is each read? | Per-read mean-Q histogram (1-Q bins, Q0-93). One line per MAPQ class. | all | Per Sequence Quality Scores |
| 3 | Where along the read does quality drop? | Mean Q (error-probability averaged) per position bin. One line per MAPQ class. | all | Per Base Sequence Quality |
| 4 | Is base composition even along the read? | % A/C/G/T/N per position bin, one panel, all reads | all | Per Base Sequence Content + Per Base N Content |
| 5 | Where along the read do errors occur? | Three panels sharing one x-axis: mismatch, insertion and deletion **rate per aligned base**. One line per mapped MAPQ class. Full page width. | aligned | none (SAMStat's differentiator) |
| 6 | How accurate are the aligned reads? | Per-read gap-compressed identity histogram, with the median marked. One line per mapped MAPQ class. | aligned | none |

**Encoding rules (one sentence each):**
- **MAPQ:** a MAPQ class is always a line with a fixed colour and dash pattern.
- **Mates:** on paired input, R1 and R2 are side-by-side panels with identical axes, never dashed lines. This removes the clash where presentation used dashes for **both** MAPQ and mate.
- **Line count:** no panel has more than 5 lines.

**FastQC check:** plots 1-4 are the FastQC modules a researcher already knows, under the same names in question form. Plots 5-6 and the MAPQ table are the only new things to learn, and they are SAMStat's reason to exist. GC, duplication, overrepresented sequences, adapters and per-tile plots are absent by design (standard-plots §6).

**Help text:** each plot has a question as its title, one always-visible sentence on what normal looks like (shape only until calibrated), and one collapsed "How to read this" block. All of this text lives in `src/report/report_help.c`.

---

## 3. Cross-topic contradictions and how they are resolved

| # | Conflict | Topics | Resolution |
|---|---|---|---|
| C1 | Plot order: all-input plots first (standard-plots) vs differentiator first and composition last (presentation) | standard-plots, presentation | **standard-plots order.** FASTA/FASTQ reports are then a clean prefix (plots 1-4) with no "not applicable" gaps in the numbering. The differentiator gets the most space, not the first slot. |
| C2 | MAPQ as a stacked-bar plot vs an HTML table | standard-plots, presentation | **HTML table, which doubles as the legend.** That is one plot fewer. |
| C3 | Quality along the read: mean Q per MAPQ class vs median+IQR | standard-plots, presentation | **Mean Q (error-probability averaged), one line per MAPQ class.** Median+IQR per class means 4 bands, which is unreadable. Median+IQR also needs a per-position quality histogram, and so the Q-cap question. |
| C4 | Quality cap: Q41 / Q60 with overflow / Q93 | long-reads D1, standard-plots, presentation | **The per-position quality matrix is replaced by two integer accumulators per bin: base count and fixed-point error-probability sum from a 94-entry lookup table.** There are no Q rows, so the cap question disappears for plot 3. The per-read Q histogram spans Q0-93. Integer accumulators keep the output byte-identical for any `-t` (the long-reads determinism rule). |
| C5 | Headline quality: mean Q over bases vs median read Q | standard-plots/presentation vs long-reads | **Mean Q over all bases, error-probability averaged**, with a fixed-point integer sum. Median read Q is readable off plot 2 and is not repeated. |
| C6 | Layout switch: first chunk, only one axis collected (long-reads) vs both axes collected and the choice made at report time (standard-plots) vs at end of run (cli-ux) | long-reads, standard-plots, cli-ux, performance | **Collect both bounded axes and choose at the end from the length histogram (median >= 1 kb means relative).** This works with stdin, needs no coordination before workers start, and handles mixed files. Absolute axis = exact positions 1-1000 (long-reads), not standard-plots' geometric 250/500 scheme. Relative axis = 100 bins. JSON carries both. Render-time `clu_data` binning is deleted, because 1000 points per line is small as SVG. |
| C7 | No flags (standard-plots) vs 5-row OK/CHECK/FAIL table plus a terminal `[OK]` status (presentation) vs a modbase CHECK row | standard-plots, presentation, modbases | **No graded checks in 2.4.** Every proposed threshold is marked UNVERIFIED in its own topic, and an uncalibrated light is the FastQC failure the maintainer wants to avoid. Only factual notes are kept (section 4). Checks can come back in 2.5 as one function, after calibration on real Dorado, Revio and NovaSeq data. |
| C8 | Base composition: one panel per MAPQ class (both) and R1 only (standard-plots D6) | standard-plots, presentation | **One panel over all reads (per mate on paired input).** Base composition describes the library, not the mapping, and the MAPQ story is told by plot 5. Four MAPQ panels times two mates would give 8 panels. *Maintainer may overrule* (paper headline example; see D-1 below). |
| C9 | Length plot traces: reads+bases (presentation) vs one line per MAPQ class (standard-plots) | standard-plots, presentation | **Reads per MAPQ class.** Bases per bin is still collected for N50 and kept in the JSON, but not drawn. |
| C10 | Missing MD hides plots 4-5 (presentation) vs mismatch from MD or `=`/X CIGAR, identity from NM (long-reads) | presentation, long-reads, cli-ux | **long-reads semantics.** The mismatch panel needs MD or X ops. The insertion and deletion panels always render. Identity needs NM or `=`/X. Records without the needed tags are counted, and the count is shown in the caption. A value is never silently zero. |
| C11 | Summary strip: `supplementary per 100` (long-reads) vs `% MAPQ >= 30` (standard-plots, presentation) | long-reads | **standard-plots' 8.** Secondary and supplementary go in a footnote as raw counts. No per-100 metric is added. |
| C12 | Terminal line: status + headline numbers (presentation) vs `wrote ... (format, reads, layout, seconds)` (cli-ux) | presentation, cli-ux | **One merged line, no status:** `samstat: sample.bam -> sample.bam.samstat.html (BAM, 12.4M reads, N50 18.2 kb, mean Q 17.1, 94.2% mapped, long layout, 41 s)`. Add `, partial: first N reads` under `-p`. |
| C13 | Default `-t`: min(4, ncores) (cli-ux) vs fixed 4, no detection (performance) | cli-ux, performance | **Fixed 4, no core detection.** `--threads` is added as the long name and `--nthreads` stays as an alias. |
| C14 | `-o/--outdir` as a new primary name (cli-ux) vs keep `-d` | cli-ux, ecosystem | **Keep `-d/--dir` as the documented name and add `mkdir -p`.** A renamed option means doc and nf-core churn for no new capability. |
| C15 | CRAM reference: `--reference` (performance) vs env vars only (cli-ux, release-cleanup) | performance, cli-ux, release-cleanup | **Env vars only (`REF_PATH`/`REF_CACHE`).** The error message names them. Add the option only after a user asks. This keeps "zero new options". |
| C16 | Deprecated options removed in 2.5 (release-cleanup) vs 3.0 (cli-ux) | release-cleanup, cli-ux | **Warn in 2.4, remove in 3.0** (removing an option is breaking). Use one warning text for all four: `samstat: warning: <opt> is deprecated and ignored`. |
| C17 | JSON `schema_version` minor bump (long-reads) vs `format_version` 2 (ecosystem) | long-reads, ecosystem, standard-plots | **`format_version` 2, one break.** Floats are printed `%.4f` and derived from integer counts (testing D1). |
| C18 | MultiQC additions: error profiles by MAPQ plus all summary values (standard-plots) vs 2 columns + one MAPQ >= 30 error linegraph (ecosystem) vs + median read Q column (long-reads) | standard-plots, ecosystem, long-reads | **ecosystem's set plus `samstat_mean_q`** (it is the headline value). There are no per-MAPQ plots in MultiQC. The error linegraph uses the same axis as the sample's HTML. |
| C19 | Reader batch: 10k records / 32 Mbases (long-reads) vs 4096 / 8 Mbases (performance) | long-reads, performance | **performance's 4096 / 8 Mbases**, because performance owns the reader. It is internal and has no option. |
| C20 | Lazy R2 allocation: in (performance) vs out (long-reads) | | **Out.** Memory is already bounded. |
| C21 | Perf gating: perf CTest label (performance, long-reads) vs manual script, no gates (testing) | performance, long-reads, testing | **testing: one `tests/perf/bench.sh`, no perf label, no ratio gates.** `scripts/benchmark.sh` is not added. **One TSAN CI job is kept** (from performance), because intra-file threading is new and ASAN does not find races. |
| C22 | Ultralong size ceilings: HTML < 500 KB / < 1 MB / < 2 MB | presentation, standard-plots, testing | **JSON < 1 MB and HTML < 1 MB** for the 1 Mb-read test. |
| C23 | Bioconda now with 2.3.0 (release-cleanup) vs after 2.4.0 (ecosystem) | release-cleanup, ecosystem | **2.3.0 now.** The autobump bot then handles 2.4.0. Drop `bzip2`, `xz` and `libcurl` from the host deps now, and re-add `bzip2` and `xz` in the same commit that enables them for CRAM. |
| C24 | System-htslib minimum: raise (performance) vs keep 1.16 (release-cleanup) | performance, release-cleanup | **Keep 1.16 unless an API in use needs more.** Now that modbases is deferred (section 5), nothing needs MN or 1.18. |
| C25 | libdeflate: auto-detect (performance) vs keep disabled (release-cleanup) | | **Keep disabled in 2.4.** Builds that differ by machine are worse than an unmeasured speed gain. Revisit with `bench.sh` numbers. |
| C26 | Header and layout wording appears in three topics | long-reads, cli-ux, presentation | **One header line:** `file · format · short/long layout (median length X) · SAMStat version`, plus the `-p` partial note when it applies. |

---

## 4. Notes (the only non-plot messages in the report)

They are shown under the summary strip only when they apply. They state facts, never verdicts, and need no thresholds:

1. Partial run: "first N reads only (`-p`); not a random sample".
2. Mismatches unavailable or partial: "mismatch profile excludes N reads without MD tag or X operations".
3. Identity unavailable or partial: "identity unavailable for N of M mapped reads (no NM tag)".
4. Layout fallout, only when the file is mixed: "N% of bases lie beyond position 1000 and are not in the position plots".

Truncated or unreadable input is not a note. It is a per-file error with exit 1, as today.

---

## 5. Cut or deferred

| Item | Source topic | Decision | Why |
|---|---|---|---|
| OK/CHECK/FAIL table, terminal `[OK]` status, threshold struct, `docs/reading-the-report.md` thresholds table, check unit tests | presentation D4/D5 | **Deferred to 2.5** | Thresholds are UNVERIFIED by the topic's own admission, and false alarms are FastQC's known failure. |
| Info lines with numeric thresholds (% mapped < 50, MAPQ<10 > 20%) | standard-plots §5 | **Cut** | These are thresholds under another name. Facts only (section 4). |
| Modified-base block: two summary rows, JSON `modifications` object, 2 MultiQC columns, fixtures | modbases | **Deferred entirely** (additive later) | No user has asked for it (modbases §1), it is biology-adjacent, and it is 120 lines plus fixtures for one maintainer. It also breaks the "same cells every run" rule. |
| MAPQ stacked bar plot | standard-plots plot 5 | **Cut**, replaced by the table | Four numbers do not need a chart. |
| Base composition split by MAPQ (4 panels) | standard-plots, presentation | **Cut** to one panel | See C8. |
| Reads+bases double trace in the length plot | presentation | **Cut** | See C9. |
| Per-plot "Data" `<details>` table | presentation D3 | **Cut** | The JSON is the data channel. It doubles the HTML size and nobody reads it. |
| `docs/reading-the-report.md` | presentation | **Deferred** | It duplicates `report_help.c`, and two copies of the help text will drift. |
| `supplementary/secondary per 100` | long-reads §1.7 | **Cut** | Raw counts in the footnote are enough. |
| `-o/--outdir` rename | cli-ux | **Cut** | See C14. Keep `mkdir -p`. |
| `--reference` | performance | **Deferred** | See C15. |
| Core auto-detection, min(4, ncores) | cli-ux | **Cut** | See C13. |
| Geometric absolute bins (250 exact, then 5%, cap 500) | standard-plots §4 | **Cut** | The 1000-exact axis already exists in the code. |
| Render-time `clu_data` interval binning | standard-plots, presentation | **Deleted** | Bins come from the collector, and ≤1000 SVG points are cheap. |
| `--plotend` and the `_end` matrices | all | **Removed** (hidden no-op until 3.0) | They grow with read length, and the relative axis shows read ends. |
| `--norm` | all | **Hidden no-op until 3.0** | The axis is chosen automatically. |
| Perf CTest label, ratio gates, `scripts/benchmark.sh` | performance, long-reads | **Cut** | See C21. |
| Linux RSS CTest | long-reads §8.3 | **Cut** | `bench.sh` measures it manually (testing D2). |
| libdeflate | performance D3 | **Deferred** | See C25. |
| Lazy R2 allocation | performance P0 | **Cut** | See C20. |
| Native MultiQC module, Galaxy, conda htslib switch, Docker image, nf-core copy in the repo | ecosystem, release-cleanup | **Out** (all topics agree) | |
| nf-core module | ecosystem | **After the 2.4.0 release, not part of the 2.4 code** | It depends on bioconda, and someone has to maintain it (ecosystem D4). |
| Protein FASTA | performance, cli-ux, release-cleanup | **Dropped** (all recommend it), with a CHANGELOG line | Retiring the tld reader removes two workarounds and a temp copy. |

**Kept although they add a little:**
- **TTY-only progress counter:** about 15 lines. A 20-minute WGS run with no output looks hung.
- **TSAN CI job.**
- **MultiQC `samstat_error` linegraph:** it is the only place the differentiator reaches a cohort view.
- **CRAM input:** it comes free with content detection. It needs bz2 and lzma enabled in the bundled htslib.

---

## 6. Ownership clarifications for the implementation workflow

- **Collector, `src/collect`:** long-reads is the single owner of fields and accumulators. It applies C4, C5, C6 and C11. standard-plots §7 and performance P0 defer to it.
- **Rendering:** presentation owns the SVG primitives and `report_help.c`. It uses the plot table in section 2 of this file, not its own D2 table.
- **Reader:** performance owns it, including batch size, content detection, stdin, CRAM and the `-p` cutoff semantics. cli-ux owns only option parsing, path checks and the terminal line.
- **JSON v2 field names:** they are reserved first (ecosystem §3), then frozen. Testing's wave 0 writes the disabled tests against them.
- **Golden regeneration:** at most twice (testing T9). The `-p` secondary/supplementary change (performance D5) and the full-length change (long-reads D5) ride in regeneration #1.

---

## 7. Decisions this review still leaves to the maintainer

Each item states the ruling this review applied. The maintainer can overturn any of them.

- **D-1 Base composition:** one panel over all reads (applied), or MAPQ small multiples as in the paper's headline figure.
- **D-2 No graded checks in 2.4:** applied. The alternative is presentation's 5-row table, with platform-dependent rows shown as info only.
- **D-3 Modbases deferred entirely:** applied. The alternative is to ship only the "inconsistent MM/ML" count as a note in section 4.
- **D-4 SVG in C vs inlined plotly-basic:** all topics recommend SVG. This review agrees, because SVG is the only option that meets the HTML < 1 MB ultralong ceiling.
- **D-5 Full read length for the length plot and N50:** applied (all topics agree).
- **D-6 Layout threshold:** median 1 kb, no override (applied).
- **D-7 Before tagging:** one real Dorado BAM and one real Revio BAM are checked by hand (testing D4). This is required before any caption text quotes a typical range.
