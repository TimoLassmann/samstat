# Topic: standard plots (final design for 2.4)

Status: final design. Base: the minimal proposal. Grafted from the researcher proposal: base composition kept, an always-visible one-line caption, the axis named in the axis title, and the error profile drawn as one figure with three panels. Items marked DECISION NEEDED are for the maintainer.

## 1. Principle

Every report contains one summary strip and **at most seven plots**, always in the same order and drawn the same way, for short or long reads, aligned or not. There is no `--long` flag and no plot options. One rule holds throughout: **on aligned input, every plot has one line per MAPQ class** (>=30, 10-29, <10, unmapped). That stratification is what makes SAMStat different (paper: "revealing differences between well and poorly mapped reads", research/paper-samstat.md). A researcher learns it once and can read every plot the same way.

Verified in code: base composition, quality and length are already stratified by MAPQ. In `src/collect/collect.c`, around lines 301 and 331-332 index with `group_size * q_idx`, and line 403 does the same for `len_dist`. The four bins are defined at `collect.c:18-21`. 2.4 keeps this structure and changes the x-axes, the units and the plot list.

## 2. Summary strip (plain HTML table, first thing on the page)

It shows the same eight values in the same order every time. It replaces the file-size/mtime header.

| Value | Definition | Aligned input only |
|---|---|---|
| Reads | primary records (secondary and supplementary counts go in a footnote) | |
| Bases | sum of read lengths | |
| N50 | from the log length histogram, which stores a base count per bin; accurate to within 1% | |
| Median length | from the same histogram | |
| Mean Q | `-10*log10(mean error probability)` over all bases | |
| % mapped | mapped primary / primary | yes |
| % MAPQ >= 30 | share of primary reads in the top MAPQ class | yes |
| Median identity | gap-compressed identity (see plot 7) | yes |

- On unaligned input the alignment columns are omitted, not shown as blank.
- Mean Q from Phred averaging is wrong for long reads. Use error-probability averaging, as nanoq does (research/nanoplot-cramino.md).

## 3. The seven plots (fixed order)

The plot title is a question. Below it is **one always-visible sentence** saying what normal looks like, with no numeric thresholds. A collapsible `<details>` block holds a longer "how to read this" text and the data table. Captions must not depend on JavaScript.

| # | Plot (title question) | Contents | Input | Why it is in the set |
|---|---|---|---|---|
| 1 | **How long are the reads?** | Read-length histogram on a log-x axis with bins about 1% wide (geometric), one line per MAPQ class, and a vertical N50 line | all | The first check anyone makes: size selection, degradation, adapter dimers, and mixtures. Unmapped reads with a bimodal length distribution was one of the paper's examples. |
| 2 | **How good is each read?** | Histogram of per-read mean Q (error-probability mean, 1-Q bins), one line per MAPQ class | all | Shows whether the run is usable. It works on binned qualities, unlike boxplots. The lower quality of unmapped reads was a paper example. |
| 3 | **Where along the read does quality drop?** | Mean Q per position bin, one line per MAPQ class | all | 3' decay, failed cycles, and where to trim. |
| 4 | **Is base composition even along the read?** | % A/C/G/T/N per position bin. Aligned input gets one small panel per MAPQ class, unaligned input one panel. No dropdown. | all | Adapters, N spikes and hexamer bias, without an adapter library. The N spike at position 40 was the paper's headline figure. |
| 5 | **How well did reads map?** | One horizontal stacked bar of the 4 MAPQ classes, with counts and % labelled directly | aligned | Context for plots 6 and 7. |
| 6 | **Why do reads map poorly? (errors along the read)** | One figure, three panels side by side: mismatch, insertion and deletion **rate per aligned base** per position bin, one line per MAPQ class | aligned | SAMStat's differentiator. Rates replace raw counts, which answers issue #12. It gets the full page width. |
| 7 | **How accurate are the aligned reads?** | Histogram of per-read gap-compressed identity, one line per MAPQ class, with the median marked | aligned | The standard long-read accuracy check. Short-read data sits near 100% and that is fine. Biostars users currently compute it by hand from NM. |

Paired input:
- R1 and R2 appear as solid and dashed lines in the **same** plot. There are never separate R1 and R2 plots.
- Plot 4 is the exception: it uses R1 only, plus a one-line note that R2 is in the JSON. This keeps the plot to four lines per panel. (DECISION NEEDED D6)

Definitions:
- **Error rate (plot 6):**
  - Per position bin, the rate is events divided by aligned bases covering that bin, in the read's orientation (`fix_orientation` already does this).
  - Mismatches need the MD tag. Code check: `get_aln_errors` reads MD, and the current caption says no mismatch profile is shown without MD (`stat_report.c:219`).
  - Without MD, the mismatch panel shows the line "mismatches not available: no MD tag". The other two panels still render.
- **Identity (plot 7):**
  - Formula: `1 - (X + I_events + D_events) / (M_matched + X + I_events + D_events)`, where X = NM - (inserted bases + deleted bases).
  - X, I and D come from NM plus CIGAR (minimap2 `de` definition: https://lh3.github.io/minimap2/minimap2.html).
  - NM is already read (`htsglue.c:165-168`) but unused, and `a->error` is uninitialised when NM is absent. That must be fixed.
  - Without NM, the plot is replaced by the line "identity not available: no NM tag".
  - Behaviour on real Dorado and Revio BAMs is UNVERIFIED. The checklist is in research/longread-metrics.md.
- **Q range:**
  - Raise `QUAL_MAX_BIN` from 41 (`collect.c:43`) to 93, the SAM maximum: QUAL is ASCII 33-126 (https://samtools.github.io/hts-specs/SAMv1.pdf).
  - The per-position quality matrix stores mean-Q accumulators (sum of error probability and count), not a 94-row histogram. Memory stays bounded.

## 4. One position-axis rule (plots 3, 4 and 6)

- The collector fills **both** bounded axes for every read, so it needs one pass and no prior knowledge of read lengths:
  - A: absolute position, exact for 1-250 bp, then about 5% geometric bins, capped at 500 bins.
  - B: relative position, 100 bins (% of read length; this is the existing `--norm` idea).
- At report time:
  - Median read length below 1 kb shows A, labelled "position (bp)".
  - Otherwise it shows B, labelled "position (% of read)".
  - The JSON carries both, so MultiQC and the user can still use the other axis.
- Memory no longer depends on the longest read. The 1 Mb read case drops from 3.8 GB RSS and a 96 MB JSON (research/code-io.md, research/code-report.md) to a fixed few MB.
- Binning happens in the collector, not at render time: `clu_data` in `plot_interval.c` and its `round(mean)` behaviour go away.
- The rate and percentage plots stay meaningful in sparse tail bins because they divide by coverage.

## 5. Flags: none in 2.4

- Plots and summary values carry no pass/warn/fail markers. FastQC's thresholds produce false FAILs on RNA-seq, amplicon and long-read data (research/fastqc.md). Researchers stop trusting a red light that is often wrong.
- The single exception is an "info" line under the summary strip, shown only when a value is unambiguous:
  - % mapped < 50
  - % MAPQ < 10 > 20
  - mismatch panel or identity plot unavailable (no MD / NM)
- These thresholds are a proposal and have not been validated on real data.
- An OK/CHECK table would give MultiQC a status column. A later release can add it as one function, if users ask for it. (DECISION NEEDED D1)

## 6. Explicitly out

- **Not added:** GC plot, per-base N plot, duplication, adapter, k-mer, overrepresented sequences, per-tile, 2D length-by-quality or length-by-identity plots, KDE/hex plots, time and channel plots, coverage, insert size, homopolymer analysis, configurable MAPQ bins, multi-sample views (left to MultiQC), and interactive toggles or dropdowns.
- **Options removed:**
  - The "at end of reads" plots and `--plotend`.
  - `--norm`, which the automatic axis replaces. Both stay accepted as no-ops with a warning for one release (research/code-cli.md).
- **Out of this topic:** the modified-base block. If one is added, it goes below plot 7 and appears only when MM/ML tags are present (research/modbases.md).

## 7. Affected files

- **`src/collect/collect.c`:**
  - Collect directly into bounded axes A and B, the log length histogram (reads and bases per bin), the per-read mean-Q histogram, the identity histogram, and mean-Q accumulators.
  - Drop the `_end` and `len_norm` matrices, and the `max_len` resizing at lines 92-116.
  - This is shared with the collector and performance topics, which own it.
- **`src/htsinterface/htsglue.c`:** set `a->error` to -1 when NM is absent; collect the CIGAR I/D event and base counts.
- **`src/plot/*` and `src/report/stat_report.c`:**
  - Build the summary strip, the 7-plot list in fixed order, captions, and the not-available lines.
  - The rendering technology (inline SVG in C or inlined Plotly) belongs to the presentation topic. (DECISION NEEDED D5)
- **`src/report/json_report.c`, `mqc_report.c`, `docs/json-format.md`, `docs/samstat.schema.json`:**
  - The same series with fixed-size arrays and a schema version bump.
  - The MultiQC output gains error profiles by MAPQ and general-stats columns for the summary-strip values.
- **`src/param/param.c`:** deprecate `--plotend` and `--norm`.

## 8. Effort (estimate, single developer)

| Piece | Days |
|---|---|
| Collector, bounded axes and new histograms | 3-4 |
| Plot list, merged error figure, summary strip, captions | 3-4 (rendering technology not included) |
| JSON, schema, MultiQC | 1-2 |
| Tests and golden regeneration | 2 |
| **Total** | **about 10-12** |

## 9. Risks

- **Golden churn.** Every report changes. Regenerate the goldens in one dedicated commit after the collector lands.
- **Changed meaning for 2.3 JSON consumers.** Counts become rates, and the length histogram moves to log bins. Mitigation: schema version bump plus a changelog entry.
- **Caption wording.** The "normal looks like" text is UNVERIFIED against real long-read data. Keep it qualitative until checked on a handful of real BAMs.
- **Binned qualities.** On NovaSeq X and similar data, the per-read Q histogram has few occupied bins. This is acceptable and needs no special mode, but the caption should mention it.
- **Mixed-length files.** Some files are mixed, for example short adapters plus long reads. The median rule picks one axis, and the other axis is still available in the JSON.

## 10. Tests

- **Ground truth from generated fixtures** (`tests/data/longread/`, research/code-tests.md):
  - N50 and median within 1%.
  - Identity from known CIGAR and NM, including an indel.
  - Mean Q computed over error probabilities (a constructed case where the Phred mean differs).
  - Q > 41 kept.
- **Fixed set:** the plot IDs and titles appear in the specified order.
  - Aligned input gives plots 1-7.
  - FASTQ input gives plots 1-4 only.
  - No MD gives the mismatch "not available" line.
  - No NM gives the identity "not available" line.
- **Axis rule:** a 150 bp median selects "position (bp)" and a 10 kb median selects "position (% of read)". Both axes are always present in the JSON.
- **Size bound:** one 1 Mb read gives JSON and HTML under 1 MB each and an RSS ceiling. The 4 MB golden skip is lifted.
- **Error rates:** a tiny SAM with known events gives exact per-bin rates for each MAPQ class.
- **No flags:** no PASS/WARN/FAIL strings appear. The info line appears exactly at its boundary conditions.
- **Unchanged:** exit codes 0/1/2.

## DECISION NEEDED

- **D1. Flags.** Recommended: none in 2.4, only info lines. Alternative: an OK/CHECK/INFO table, which would give MultiQC a status column.
- **D2. Base composition (plot 4).** Recommended: keep it, because it was the paper's headline example and is already collected. Alternative: cut it to reach 6 plots.
- **D3. Length definition (plot 1 and N50).** Recommended: full read length. Today the code uses clip-trimmed length (`collect.c:403`).
- **D4. Long-read switch.** The proposed threshold is a median length of 1 kb, with 100 relative bins. Should there be no override at all (recommended)?
- **D5. Rendering.** Hand-written inline SVG (recommended by research/presentation.md) or inlined plotly-basic (about 1.2 MB, per https://github.com/plotly/plotly.js/blob/master/dist/README.md). This is decided in the presentation topic.
- **D6. Paired input.** Should plot 4 show only R1 (recommended, R2 stays in the JSON), or R1 and R2 panels?
