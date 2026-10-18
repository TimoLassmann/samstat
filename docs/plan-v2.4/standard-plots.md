# SAMStat 2.4: the standard report (exact specification)

Status: binding spec for the implementation workflow. Sources: `simplicity-review.md`, `topics/standard-plots.md`, `topics/long-reads.md`, `topics/presentation.md`, `topics/ecosystem.md`. Precedence: **for exact text, data definitions and rendering rules this file wins; for scope, budget and cuts `simplicity-review.md` wins.** Where this file adds detail the topics left open, the choice is listed in `decisions.md` section F. Changes made after the plan review are logged in `review-log.md`.

The rule: a FastQC user reads the whole report in 30 seconds the first time. Same page, same order, same encodings for every input. No options change the report.

---

## 1. Page order (fixed)

1. **Header line** (one line).
2. **Summary strip** (8 values on aligned input, 4 on unaligned input). A file with no reads shows neither the strip nor anything below the header except the banner `This file contains no reads.` and the footer.
3. **Notes** (factual lines, only when they apply).
4. **MAPQ table** (aligned input only). It is also the legend for every MAPQ line on the page.
5. **Plots.** Unaligned input shows plots 1-4 and nothing after them (no N/A lines for 5-6). Aligned input shows plots 1-6. A plot that applies to the input type but has no usable data (FASTA without qualities for plots 2-3; no mapped reads, no MD/X or no NM for 5-6) is replaced in place by a one-line N/A note; numbering never changes.
6. **Footer:** SAMStat version, input basename, citation (Lassmann T., "SAMStat 2: quality control for next generation sequencing data", Bioinformatics 39(1), btad019, 2023, doi:10.1093/bioinformatics/btad019, as in `CITATION.cff`), raw secondary and supplementary counts. **No command line** in the footer: it would leak absolute paths and break byte-identity across `-t` (`decisions.md` F10).

There is no flag, status light, dropdown, tab, zoom, hover tooltip or JavaScript. Output: one self-contained HTML file `<input>.samstat.html` with inline SVG (`decisions.md` A9). The file makes zero network requests.

**HTML size ceiling:** < 1 MB for both layouts, tested on the 1 Mb single read (long layout) and on a paired 999 bp aligned file (short layout, the worst case for point count).

### 1.1 Header line

`<file> · <format> · <short|long> reads (median <X>) · SAMStat <version>`

- When no read has a sequence (`total_bases` = 0), the reads clause is left out (the layout is a guess from zero lengths). The word "layout" and the 1 kb switch are not in the header: the switch is explained in the collapsed note of plot 1 ("Files whose median read length is 1 kb or more are drawn in the long-read style: ...").
- `<format>` is the htslib-detected format: `SAM`, `BAM`, `CRAM`, `FASTQ`, `FASTA` (with `gzip` where it applies).
- Under `-p N`, append: ` · partial: first N reads only (-p); not a random sample`. 2.3 already renders a partial warning section (`src/report/stat_report.c:58-59`, `samstat_partial_warning_section`); 2.4 moves it into the header line.

### 1.2 Summary strip (same order every time)

| # | Label | Definition | JSON key (`summary.`) | Aligned only |
|---|---|---|---|---|
| 1 | Reads | Primary records. Paired input labels it `Reads (R1+R2)`. Secondary and supplementary records are excluded from every statistic. | `total_reads` | |
| 2 | Bases | Sum of full read length (`l_qseq`, soft-clipped bases included); shown in bases (`9,876 bp`, `89.5 kb`, `1.0 Mb`) with the exact count as the small line when the display is rounded | `total_bases` | |
| 3 | N50 | Smallest length L such that reads of length >= L hold >= 50% of bases. From the length histogram (section 3.1). Small line: `half of all bases are in reads at least this long`. | `n50` | |
| 4 | Mean Q | Error-probability average over **all bases with a quality value, clipped bases included** (section 3.3). Omitted for FASTA. Small line: `average base quality (Q30 = 1 error in 1000 bases)` | `mean_q` | |
| 5 | % mapped | Mapped primary / primary x 100 | `pct_mapped` | yes |
| 6 | % MAPQ >= 30 | Primary reads with MAPQ >= 30 / primary x 100. Small line: `reads placed with high confidence` | `pct_mapq_ge30` | yes |
| 7 | % bases clipped | Soft-clipped bases of mapped primary reads / all bases of mapped primary reads x 100. Small line: `ends the aligner trimmed off` | `pct_clipped` | yes |
| 8 | Median identity | Median per-read gap-compressed identity (plot 6). Small line: `how closely reads match the reference` | `median_identity` | yes |

- Median length is not a strip value: it is in the header line and in JSON `summary.median_length`.
- Unaligned input: values 5-8 are omitted, not shown as 0 or blank.
- No sequence in any read (`total_bases` = 0, e.g. SEQ is `*`): Bases and N50 read `n/a`, never `0`; the note in 1.3 says why.
- The three percentage tiles (5, 6, 7) and Median identity carry a `%` sign.
- Number formatting in HTML: counts with thousands separators and SI suffix (`12.4 M`), lengths in bp below 10 kb, in kb up to 1 Mb (`18.2 kb`) and in Mb above (`1.0 Mb`, the unit of the Bases tile; the same in the N50 label of plot 1), Q to one decimal, percentages to one decimal.
- Length values (N50, median length) are exact below 1000 bp and within 1% above; above 1 kb the HTML prefixes `~` unless the value is known exactly (section 3.1: one read in the bin, or every read has the same length).

### 1.3 Notes (facts, never verdicts, no thresholds)

Shown under the strip only when the condition holds. This table is the only source for note wording and conditions (it supersedes `simplicity-review.md` §4). Exact text:

| Condition | Text |
|---|---|
| `-p` used | `Partial run: first N reads only (-p); not a random sample.` |
| no read has a sequence | `Reads have no stored sequence (SEQ is *): length, quality and base composition are not available.` |
| 50 to 99 reads (and some sequence) | `Few reads (N): per-position plots are noisy, and positions covered by fewer than 50 reads are not drawn.` |
| fewer than 50 reads (and some sequence) | `Only N reads: per-position plots (3-5) are not reliable and are hidden.` (`read` for N = 1) |
| `reads_without_md > 0` | `Mismatch profile excludes N of M mapped reads (no MD tag and no =/X operations in CIGAR). Add MD with: samtools calmd -b in.bam ref.fa > out.bam (minimap2: use --MD).` |
| `identity_na > 0` | `Identity unavailable for N of M mapped reads (no NM tag and no =/X operations).` |
| short layout and `bases_beyond_position_axis > 0` | `P% of bases lie beyond position 1000 and are not in plots 3-5.` |

`minimap2 --MD` per https://lh3.github.io/minimap2/minimap2.html (UNVERIFIED, page not fetched). Truncated or unreadable input is not a note; it is a per-file error, exit 1, no report.

### 1.4 MAPQ table "How well did reads map?" (aligned input only)

Caption (one sentence): `MAPQ is the aligner's confidence that a read is placed correctly; aligners scale it differently.`

Orange is darkened from Okabe-Ito `#E69F00` (2.25:1 on white) to `#B86E00` (4.0:1): every series colour has at least 3:1 contrast on white (WCAG 1.4.11; tested). The class name appears only in the row header; the line sample's accessible name describes its style (`dashed orange line`).

Four rows in fixed order. Columns: line sample (colour + dash, drawn as a short inline SVG), class label, reads, % of primary reads. Paired input adds R1 and R2 count columns.

| id (JSON, frozen) | Label (2.4) | Colour (Okabe-Ito) | Dash |
|---|---|---|---|
| `mapq_ge30` | MAPQ >= 30 | blue `#0072B2` | solid |
| `mapq_10_29` | MAPQ 10-29 | dark orange `#B86E00` | dashed `8 4` |
| `mapq_lt10` | MAPQ < 10 | vermillion `#D55E00` | dotted `3 4` |
| `unmapped` | Unmapped | grey `#777777` | dash-dot `12 4 2 4` |

- MAPQ 255 ("unavailable" in the SAM spec, https://samtools.github.io/hts-specs/SAMv1.pdf, UNVERIFIED; STAR uses it for unique hits) is classed as `>= 30`, unchanged from 2.3 (`collect.c` `get_mapqual_bins`: every value >= 30).
- Okabe-Ito hex values: https://jfly.uni-koeln.de/color/ (UNVERIFIED, page not fetched; values as widely reproduced). Grey is not an Okabe-Ito colour; `#777777` is chosen for >= 4.5:1 contrast on white (UNVERIFIED, check with a contrast tool during WP4.2b).

Unaligned input: the table is absent and plots 1-4 draw one line, `all reads`, solid black `#000000`.

---

## 2. Encoding rules (apply to every plot)

1. **MAPQ class = one line** with the fixed colour and dash above. On aligned input plots 1-3, 5 (insertion and deletion rows) and 6 have one line per class that has data (unmapped is not drawn in plots 5 and 6). Plot 4 and the mismatch rows of plot 5 are the exception: there colour means the base, so each MAPQ class gets its own row instead of its own line (rule 9).
2. **Mates:** paired input draws R1 and R2 as side-by-side panels with identical x and y axes and the titles `Read 1` / `Read 2`. Never extra lines.
3. **Max 5 lines per panel.**
4. **No direct line labels.** MAPQ lines are identified only by the MAPQ table. Plots 4 and 5 each have one fixed inline legend row (A C G T N swatches) in their caption area; plot 4 and the plot 5 mismatch rows are labelled by MAPQ class in the row title.
5. Each plot: title (a question), one always-visible caption sentence (above the figure), one always-visible `What good looks like:` line (under the figure; it describes the ideal, not a verdict on this file), one collapsed `<details>` `How to read plot N` with `Worth a look:`, `Often benign:` and, for plots 2 and 6, `How it is calculated:`. All text lives in `src/report/report_help.c` (section 5). The caption of plot 3 names MAPQ classes only for aligned input.
6. Every `<svg>` has `role="img"`, `<title>`, `<desc>`. Each series is one `<path>`/`<polyline>` with one series-level `<title>`; **no per-point elements or tooltips**. The JSON is the data channel. Lines are 2 units wide. A histogram bin with no drawn neighbour is drawn at least 6 units wide and 1 unit high, so a single read on a log axis is still visible. A log axis has at least three labelled ticks (200 k, 500 k, 1 M for a lone bin), a count axis is at least 0-3 and uses whole-number ticks.
7. Figures are 400 units wide with 15 unit text and scale to their column (a figure up to 32 rem wide, a Read 1 / Read 2 pair side by side, one column below 40 rem), so text is >= 12 px at 360 px width. System font stack, >= 12 px text, light theme only, layout works at 360 px width, one `@media print` block (forces `<details>` open, `break-inside: avoid`, A4 portrait).
8. The base palette (plots 4 and 5) is the IGV nucleotide palette and shares no hex value with the MAPQ palette. Some hues are close (C blue / MAPQ >= 30 blue, G orange / MAPQ 10-29 and < 10); they never meet in one row: bases are filled areas or base-coloured lines in rows titled by MAPQ class, MAPQ classes are lines in the plot 1-3, 6 and indel panels.
9. **Base x MAPQ = rows.** Where a plot shows bases (plot 4, plot 5 mismatch rows), each MAPQ class that has data is one short row, in the MAPQ table order, titled with the class label (`MAPQ >= 30`, ..., `Unmapped`); the five bases are the colours inside the row. A class with no reads in either mate has no row. Rows share the x-axis; paired input shows R1 and R2 as columns, as everywhere else. Unaligned input has one row, `all reads`. This restores the contrast both SAMStat papers are built on (decisions A3, F3 revised 2026-10-18).

---

## 3. Data model (collector, `src/collect/collect.{c,h}`)

All accumulators are fixed-size integer arrays, independent of read length. Every statistic is primary-only (secondary/supplementary are skipped in `read_bam_chunk`, `src/htsinterface/htsglue.c:82-87`). Indices: `m` = MAPQ class (4), `r` = mate (2).

Layout in memory: per position bin, one struct holding all accumulators of that bin for one `m`, `r` (array of structs), so a base touches one or two cache lines.

### 3.0 Input contract (reader to collector)

- **SAM/BAM/CRAM:** records carry `aln_data`. `has_sq` = the header has at least one `@SQ` line (read once on the main thread before workers start).
- **FASTA/FASTQ** (`hts_get_format(fp)->format` is `fasta_format` or `fastq_format`): records are delivered **without** `aln_data` (htslib marks them `BAM_FUNMAP` and sets `BAM_FPAIRED|BAM_FREAD1|BAM_FREAD2` from `/1`, `/2` name suffixes; the reader ignores those flags). Every read is R1, as in 2.3. `parse_alignment` is not called.
- **Empty input:** a zero-byte or whitespace-only input (after decompression) is a valid empty input: exit 0, one warning, `total_reads = 0`. Checked before format detection. A FASTA/FASTQ record with an empty sequence is not counted (keeps the 2.3 `fastx` goldens). Only non-empty input whose content htslib cannot recognise fails that input with exit 1.
- **Per-record decisions, never per chunk:** a record contributes to quality accumulators iff its QUAL is present (not 0xFF); to alignment accumulators iff it has `aln_data` and is mapped. There is no `collect_type` chosen from the first chunk.

### 3.1 Length histogram (plot 1, N50, median)

- Length = full `l_qseq` (clipped bases included). Replaces the 2.3 clip-trimmed length (`collect.c:403`).
- **One master grid.** `g(e0)` = the integer sequence `e[0]=e0; e[j+1] = max(e[j]+1, (e[j]*101 + 99)/100)` up to the first edge >= 10,000,000.
  - **Storage edges** (`length_edges`, JSON): exact lengths 0..999 (bin i = length i), then `g(1000)`, plus one overflow bin. 1922 bins. Lookup by binary search. No floating point.
  - **Display edges** (plot 1 long layout, MultiQC): `0`, the edges of `g(1)` below 1000 (1..100 exact, then 101, 103, ..., 986, 996), then exactly the storage edges from 1000 upwards. 1216 bins. Every display bin is a union of whole storage bins, so summing is exact (no sawtooth).
  - Plot 1 in short layout uses the storage bins directly (1 bp below 1000).
- Per bin, per `m`, per `r`: `reads` and `bases` (sum of lengths).
- Exact scalars kept as well: `total_bases`, `min_length` (initialised to `UINT64_MAX`; reported as 0 when no reads), `max_length`.
- `n50` and `median_length` (JSON) are reported as the **lower edge** of the storage bin where the cumulative fraction (bases for N50, reads for median, counted from the longest/shortest end respectively) first reaches 50%. Exact below 1000 bp.
- **Display values** (`n50_best`, `median_best`, with `*_exact` flags; used by HTML, MultiQC and the stderr line, not by JSON): equal to the lower edge, except that a bin holding exactly one read gives that read's length (exact) and an edge is never shown below `min_length` or above `max_length`; when `min_length` = `max_length` the value is exact. A single 1 Mb read therefore reads `1,000.0 kb`, not `~995.6 kb`. JSON keeps the bin edges until the next golden regeneration (open item for the maintainer).

### 3.2 Position axes (plots 3, 4, 5)

Both axes are filled for every read in one pass; the report picks one at the end of the run. On aligned input the position axes cover **the aligned part of the read**: clipped bases are excluded, and the axis label says so.

| Axis | Bins | Mapping of read position i (0-based, read orientation via `fix_orientation`, clipped bases excluded), aligned length L | Label |
|---|---|---|---|
| absolute | 1000 | bin = i for i < 1000; otherwise not counted, added to `bases_beyond_position_axis` | `position in read (bp)` |
| relative | 100 | bin = (100*i)/L, computed by an incremental walk over bin boundaries per read (no division per base) | `position in read (% of length)` |

On aligned input the plot labels append `, aligned part`.

**Layout rule:** `layout = "long"` if `median_length >= 1000`, else `"short"`. Decided once, after all records are read (works for stdin). No override. Short layout shows the absolute axis, long layout the relative axis. JSON carries both.

**x range:** short layout x = 1 to the last bin backed by at least 50 bases (over all classes, per mate; R1 and R2 panels share the range). Long layout x = 0-100%.

**Thin positions:** in plots 3-5 a bin backed by fewer than 50 bases (`REPORT_MIN_POS_BASES`; all classes pooled, per mate; on the absolute axis that is fewer than 50 reads) is a gap. **Too few reads:** an input with fewer than 50 reads has no per-position view at all: plots 3, 4 and 5 are each replaced by the N/A line `Too few reads (N) for a per-position plot (needs at least 50).` (plots 1, 2 and 6 are drawn). The stated rule is the rule that is applied.

`bases_beyond_position_axis` (JSON and the note in 1.3) counts bases beyond the **absolute** axis only; in the long layout (relative axis) no base is lost and the field is not shown as a note.

Per bin, per `m`, per `r`:

| Accumulator | Meaning |
|---|---|
| `composition.{A,C,G,T,N}` | read bases; anything not ACGT counts as N |
| `qual_bases` | bases with a quality value (QUAL not `*`/0xFF) |
| `qual_error_sum` | sum of `E[q]` over those bases (3.3) |
| `aligned` | aligned read bases (CIGAR M, =, X, I) of mapped reads; denominator for insertions and deletions |
| `aligned_md` | same, only for reads that have an MD tag or X operations; denominator for mismatches |
| `mismatches` | mismatching bases (from MD, or CIGAR X when MD is absent) |
| `mismatch_base.{A,C,G,T,N}` | the same mismatches split by the **read** base at that position (read orientation, after `fix_orientation`); anything not ACGT counts as N. Their sum equals `mismatches` in every bin. The read base, not the reference base, so it works for CIGAR X without MD and shows N calls and bases added by library preparation (2011 paper: 5' added G in CAGE) |
| `insertions` | insertion events starting at this read position (one per gap opening) |
| `deletions` | deletion events at this read position (one per gap opening, 2.3 semantics) |

Memory: 17 rows x 1100 bins x 4 x 2 x 8 B = about 1.2 MB per collection (was 12 rows, 0.85 MB, before `mismatch_base`).

### 3.3 Quality, in integer fixed point

- Lookup table `E[q] = round(10^(-q/10) * 2^30)` for q = 0..93, generated at build time (like `convert_tables.h` from `src/convert_tables.c`) so no `pow` runs at collect time. `E[0] = 1073741824`.
- Phred values above 93 are clamped to 93 (SAM maximum: QUAL is ASCII 33-126, https://samtools.github.io/hts-specs/SAMv1.pdf, via topic research, not re-checked: UNVERIFIED). This replaces the Q41 cap (`QUAL_MAX_BIN`, `collect.c:43`).
- **Scalars over all bases (clipped included):** `qual_bases_all` (`uint64_t`) and `qual_error_sum_all` (`unsigned __int128`, internal only, not written to JSON). The 128-bit sum removes the overflow risk (a 1 Tbase run at Q15 sums to about 3.4e19, above `UINT64_MAX`).
- Per-bin `qual_error_sum` stays `uint64_t`: it overflows only after about 1.7e10 Q0 bases in one bin.
- Precision: `E[q]` rounds to 1 near Q90; mean Q values above about Q80 are coarse. Accepted.
- **Per-read mean Q** (plot 2): over all bases of the read with a quality (clipped included). Per read, `s = sum E[q]`, `n = bases with quality`, `mean = s / n` (integer division). Bin = the largest q with `E[q] >= mean`, found by binary search on the table (integer comparisons, deterministic). Bin q means "read mean Q in [q, q+1)". 94 bins per `m`, `r`. Reads without qualities are not binned.
- **Summary `mean_q`:** `-10*log10(qual_error_sum_all / qual_bases_all / 2^30)` in `double`. **Per-bin mean Q** (plot 3): the same formula on the bin's sums. Both computed at report time from integer sums.

### 3.4 Identity (plot 6, median identity)

Gap-compressed identity, minimap2 `de` convention (https://lh3.github.io/minimap2/minimap2.html, via topic research: UNVERIFIED):

- From the CIGAR (one pass over operations, in `parse_alignment`, `src/sambamparse/sam_bam_parse.c`): `M` = bases in M, =, X ops; `I_ev`, `I_bases`, `D_ev`, `D_bases` = insertion and deletion events and bases; `S_bases` = soft-clipped bases (for `pct_clipped`). N, H, P ops are ignored.
- Mismatch count `X`, in order of preference:
  1. NM present: `X = NM - I_bases - D_bases`.
  2. No NM, CIGAR has = or X ops: `X` = bases in X ops.
  3. Otherwise: not available; count in `identity_na`.
  - `X < 0` (NM inconsistent with CIGAR) also counts as not available.
- Identity = `(M - X) / (M + I_ev + D_ev)`. Clipped bases do not count. Bin = `(1000 * (M - X)) / (M + I_ev + D_ev)` (integer division), 1001 bins: bin b = identity in [b/10 %, (b+1)/10 %), bin 1000 = exactly 100%.
- Per `m` (mapped classes only), per `r`. Unmapped reads are not binned.
- `median_identity` = lower edge of the median bin, in percent.
- Bug fixed in passing: `a->error` (NM, `struct aln_data` in `src/htsinterface/htsglue.h:24`) is set only when NM is present (`htsglue.c:165-169`), so a reused buffer slot inherits a stale value. Reset to -1 per record.

### 3.5 Missing MD

- 2.3 adds mismatches only from the MD walk; `parse_alignment` sets `genome[rp] = read[rp]` for M, = and X (`sam_bam_parse.c:158-166`), so a BAM without MD silently shows zero mismatches. 2.4: X ops count as mismatches when there is no MD. MD still wins when present.
- Mapped reads with neither MD nor X ops: counted in `reads_without_md`, contribute nothing to `mismatches` or `aligned_md`. Insertions and deletions come from the CIGAR and are always counted.

### 3.6 Scalars and `is_aligned`

`total_reads`, `read1_reads`, `read2_reads`, `paired`, `properly_paired`, `secondary`, `supplementary` (2.3 meaning; counted only in the reader thread, under `-p` only up to the cutoff, `decisions.md` B5), `total_bases`, `min_length`, `max_length`, `mapped_bases`, `clipped_bases`, `qual_bases_all`, `qual_error_sum_all`, `reads_without_md`, `identity_reads` (reads binned in plot 6), `identity_na`, `bases_beyond_position_axis`, `has_sq`, `is_partial_report`.

**`is_aligned`** (computed in finalise, never set by input type): `is_aligned = has_sq || (at least one primary record is mapped)`. Always 0 for FASTA/FASTQ.
- Unaligned BAM from Dorado or PacBio (no `@SQ`, all records flag 4): unaligned layout (plots 1-4, no MAPQ table).
- `@SQ` present but no read mapped: aligned layout; the MAPQ table shows 100% unmapped; plots 5 and 6 show the N/A note `No mapped reads.`

**Merge semantics** (`stat_collection_merge(dst, src)`, used by intra-file threading): every counter and accumulator is summed; `has_sq` and `is_partial_report` are OR-ed; `min_length` is min, `max_length` is max. A collection that saw zero reads is the identity element. Output is therefore byte-identical for any `-t`.

### 3.7 Finalise

One function in `collect.c` derives every summary value (`n50`, `median_length`, `mean_q`, `pct_mapped`, `pct_mapq_ge30`, `pct_clipped`, `median_identity`, `is_aligned`, `layout`). HTML, JSON, MultiQC and the stderr line read these values; nothing recomputes them.

---

## 4. The six plots

Common: plot title = question; caption and `What good looks like:` line always visible; "How to read plot N" collapsed. Exact text in section 5. Plots 1-4 apply to all inputs, 5-6 to aligned input only (absent, not N/A, on unaligned input).

### Plot 1: How long are the reads?

| | |
|---|---|
| Data | `read{1,2}.length.reads[m]`: storage bins (short layout) or display bins (long layout), section 3.1 |
| x | short layout: read length, **linear**; long layout: read length, **log scale**, ticks at 1-2-5 decades (`100`, `1 k`, `10 k`, ...). Range = min to max non-empty bin, padded to at least 10 bp (linear) or one decade (log) centred on the data; a log axis always has at least three labelled ticks (the 2 and the 5 of the decade are labelled when fewer than three decades would be) |
| y | reads per bin, linear, from 0 |
| Lines | one per MAPQ class (all reads on unaligned input) |
| Marks | vertical dashed black line at N50, labelled `N50 <value>` |
| N/A note | no sequence in any read: `Not available: the input has no sequence.` |
| FastQC equivalent | Sequence Length Distribution |

A single-length short-read library shows one honest spike in the middle of a 10 bp window. That is expected.

### Plot 2: How good is each read?

| | |
|---|---|
| Data | `read{1,2}.read_mean_q[m]` |
| x | per-read mean Q (Phred), linear, 0 to max non-empty bin + 1. Label `average base quality per read (Q), 1-Q bins`; a short unlabelled tick at every bin edge (bins are `[q, q+1)`, drawn edge to edge, not centred on a tick) when edges are at least 3 units apart |
| y | reads per Q bin |
| Lines | one per MAPQ class |
| N/A note | FASTA or no qualities: `Not available: the input has no base qualities.` |
| FastQC equivalent | Per Sequence Quality Scores (different averaging) |

### Plot 3: Where along the read does quality drop?

| | |
|---|---|
| Data | per bin `qual_error_sum / qual_bases` on the layout axis, converted to Q (3.3) |
| x | layout axis and x range (3.2) |
| y | mean Q, linear, 0 to the next multiple of 10 above the max |
| Lines | one per MAPQ class. A bin with no base in a class, or a thin bin (3.2), is a gap in the line, never 0 |
| N/A note | as plot 2 |
| FastQC equivalent | Per Base Sequence Quality (mean only, different averaging; no box plots) |

### Plot 4: Is base composition even along the read?

| | |
|---|---|
| Data | `composition` per MAPQ class, per mate, as % of the bin's bases in that class |
| x | layout axis and x range (3.2) |
| y | % of bases, 0-100 |
| Lines | IGV nucleotide colours (revised 2026-10-18; 2.4.0 used green/purple/black/brown/blue): A green `#009600` solid, C blue `#0000FF` dashed `6 3`, G orange `#D17105` dotted `2 3`, T red `#FF0000` dash-dot `8 3 2 3`, N grey `#808080` long-dashed `12 3`. Every base has its own dash pattern (WCAG 1.4.1), so telling green A from red T never rests on colour alone in plot 4; every colour has at least 3:1 contrast on white (WCAG 1.4.11; G is the lowest at 3.5:1). None of the hex values is in the MAPQ palette. |
| Legend | one inline row of five swatches A C G T N above the panels |
| Panels | one row per MAPQ class with reads (aligned input: up to 4, `Unmapped` included, because unmapped and poorly mapped reads are where N calls, adapters and added bases show), each row titled with the class label; R1, R2 as columns. Unaligned input: one row, `all reads`. Rows are short (about half the height of a plot 3 panel) and share the x-axis (encoding rule 9, `decisions.md` A3 revised) |
| Thin bins | a bin covered by fewer than 50 reads **of that class** is a gap in that row (the rule of plots 3 and 5; on the absolute axis that is 50 bases, on the relative axis 50 reads, so a few long reads do not fill a row) |
| N/A note | no sequence in any read: `Not available: the input has no sequence.` |
| FastQC equivalent | Per Base Sequence Content + Per Base N Content |

### Plot 5: Where along the read do errors occur?

SAMStat's differentiator. Full page width, the tallest figure.

| | |
|---|---|
| Data | per bin: `mismatch_base.{A,C,G,T,N} / aligned_md` per mapped class; `insertions / aligned`, `deletions / aligned` |
| Layout | rows stacked vertically, sharing one x-axis: one **mismatch row per mapped MAPQ class with data** (`Mismatches, MAPQ >= 30`, `Mismatches, MAPQ 10-29`, `Mismatches, MAPQ < 10`), then Insertions, then Deletions. Paired: two columns (R1, R2) |
| Mismatch rows | stacked areas of the five per-base rates, bottom to top A, C, G, T, N, in the plot 4 base colours (fill, no dash); the top edge of the stack is the total mismatch rate of that class. A bin is a gap where fewer than 50 reads of that class cover it (the rule of the lines; on the relative axis a few long reads would otherwise fill every bin) or it has fewer than 50 `aligned_md` bases in that class. Legend: the A C G T N swatch row (encoding rule 4) |
| x | layout axis and x range (3.2) |
| y | % of aligned bases, linear from 0. Each row has its own y-scale (each mismatch row too: the point is where along the read each class goes wrong, and a low-MAPQ row can be 10x the MAPQ >= 30 row), **shared between the R1 and R2 columns**. A row that is zero at every position keeps a whole-number scale (0, 2, 4) and is labelled `none detected` inside the panel |
| Lines | insertion and deletion rows: one per mapped MAPQ class (>= 30, 10-29, < 10). Bins with `aligned` = 0 are gaps |
| Partial tags | the mismatch row caption appends the "excludes N of M reads" note when `reads_without_md > 0` |
| N/A note | no mapped read: `No mapped reads.` Fewer than 50 reads: `Too few reads (N) for a per-position plot (needs at least 50).` No read has MD or X: no mismatch rows; in their place `Mismatches not available: no MD tag. Add one with samtools calmd -b in.bam ref.fa > out.bam (minimap2: --MD).`; insertion and deletion rows still render |

Rates replace the 2.3 raw counts (GitHub issue #12).

### Plot 6: How accurate are the aligned reads?

| | |
|---|---|
| Data | `read{1,2}.identity[m]`, mapped classes |
| x | identity %, linear, from `floor` of the 1st-percentile bin (whole %, at most 95: the axis is at least 5 points wide, whole-number ticks) to 100; no tick, bar or axis beyond 100 . Label `identity to the reference (%), <w>-point bins`, w = the bin width; a short unlabelled tick at every bin edge when edges are at least 3 units apart |
| y | reads per bin: 0.1% bins in long layout; 1% bins in short layout (10 stored bins summed); the stored exact-100% bin joins the last bin, which ends at 100 (`[99, 100]` or `[99.9, 100]`), because a 150 bp read can only take identity steps of about 0.7% |
| Lines | one per mapped MAPQ class |
| Marks | vertical dashed black line at the median, labelled `median <value>%` |
| N/A note | no mapped read: `No mapped reads.` No usable read: `Identity not available: no NM tag and no =/X operations.` |

---

## 5. Report text (goes verbatim into `src/report/report_help.c`)

Captions describe shape only. No numeric "typical range" appears until checked on real Dorado, Revio and NovaSeq data (`decisions.md` E4). Each plot shows its `What good looks like:` line under the figure; the collapsed "How to read plot N" block holds `Worth a look:`, `Often benign:` and, where a value is calculated in a way a reader may not expect, `How it is calculated:`. Plain words: MAPQ is explained in the table caption and in the tile captions, clipping as "ends the aligner trimmed off", chimeras as "reads joined from two places".

**Plot 1. How long are the reads?**
Caption: `Reads per length bin (log scale for long reads). The dashed line is N50: half of all bases are in reads at least this long.`
- What good looks like: one peak where the library preparation or sequencing kit puts it.
- Worth a look: a second peak of short reads (adapter dimers, degraded input, failed reads), especially if those reads are unmapped.
- Often benign: trimmed short reads form a ramp below the full read length; amplicon libraries show one peak per amplicon.
- How it is calculated: files whose median read length is 1 kb or more are drawn in the long-read style: a log length axis here, and positions as % of read length in plots 3-5.

**Plot 2. How good is each read?**
Caption: `Reads by their average base quality (higher is better; Q20 = 1 error in 100 bases, Q30 = 1 in 1000).`
- What good looks like: one peak at high quality.
- Worth a look: a second, low-quality peak, especially if it matches the poorly mapped or unmapped reads.
- Often benign: instruments that bin qualities give a few spikes instead of a smooth curve.
- How it is calculated: base qualities are averaged as error probabilities, so a few bad bases pull a read down and the values read lower than a plain average of Q scores, which is what FastQC reports. Mean Q at the top of the page and plot 3 are averaged in the same way.

**Plot 3. Where along the read does quality drop?**
Caption, aligned input: `Mean base quality at each position along the read, one line per MAPQ class (see the table above).` Unaligned input: `Mean base quality at each position along the read.`
- What good looks like: a flat line or a slow decline towards the end of the read.
- Worth a look: a sudden drop at one position (a failed cycle) or a steep decline that suggests trimming.
- Often benign: the first few positions are slightly lower; long-read files use % of read length as the axis, which dilutes effects in the first few bases.

**Plot 4. Is base composition even along the read?**
Caption, aligned input: `Percentage of A, C, G, T and N at each position along the read, one row per MAPQ class.` Unaligned input: `Percentage of A, C, G, T and N at each position along the read.`
- What good looks like: roughly parallel, flat lines, and the same picture in every row.
- Worth a look, unaligned input (no MAPQ rows): `An N spike at one position, or a change in composition towards the end of the read (in unaligned input often adapter read-through; aligners usually cut adapters off the alignment, so in aligned input they are not shown).`
- Worth a look, aligned input: a row that differs from the MAPQ >= 30 row, such as more N or one base piling up at the start of unmapped or poorly mapped reads (the reason they did not map); an N spike at one position, or a change in composition towards the end of the read (in unaligned input often adapter read-through; aligners usually cut adapters off the alignment, so in aligned input they are not shown).
- Often benign: uneven first positions in random-primed libraries; amplicon and bisulfite libraries are uneven by design.

**Plot 5. Where along the read do errors occur?**
Caption: `Mismatches per aligned base at each position, coloured by the read base and one row per MAPQ class; then insertions and deletions, one line per MAPQ class. Each panel has its own vertical scale.`
- What good looks like: low, flat rates, with MAPQ >= 30 the lowest, and no single base dominating the mismatches.
- Worth a look: one base dominating the mismatches at a few positions (for example G at the first positions: bases added during library preparation; N: failed base calls; G towards the 3' end on two-colour instruments: no signal), especially if it reaches further into the read in the low-MAPQ rows; rates that rise towards one end (trim there); errors piled up at the ends of poorly mapped reads (adapters, reads joined from two places); high mismatch rates in well-mapped reads (wrong reference or contamination).
- Often benign: low-MAPQ reads are often ambiguous (repeats, reads that fit several places), not necessarily wrong. Long reads carry more insertions and deletions than short reads. A few very long reads can dominate a line.

**Plot 6. How accurate are the aligned reads?**
Caption: `Reads by how closely they match the reference (100% = perfect). The dashed line is the median.`
- What good looks like: a narrow peak close to 100%.
- Worth a look: a broad distribution or a second peak at low identity (contamination, wrong reference, reads joined from two places); a high % bases clipped at the top of the page with high identity here (reads joined from two places, or adapters hidden by clipping).
- Often benign: real differences between the sample and the reference lower identity slightly; short reads pile up at 100%.
- How it is calculated: identity = matching bases / aligned columns, each gap counted once; clipped bases are not counted.

## 6. stderr (one line per input, no status)

```
samstat: sample.bam -> sample.bam.samstat.html (BAM, 12.4M reads, N50 18.2 kb, mean Q 17.1, 94.2% mapped, long layout, 41 s)
```

- Unaligned input omits `% mapped`. `-p` appends `, partial: first N reads`.
- N50 is the display value of section 3.1 and is left out when no read has a sequence.
- A progress counter (every ~1 M reads) only when `isatty(STDERR_FILENO)`.
- A write failure for an output file names the remedy: `samstat: error: cannot write <path>: <strerror>; use -d DIR to write elsewhere.`

---

## 7. JSON `format_version` 2 (frozen for 2.x after 2.4.0)

Rules carried over from v1: deterministic, no timestamps or absolute paths, fixed key order, trailing zeros trimmed (an index past the end means 0), every key always present. New: floats appear only in `summary`, printed `%.4f`, computed from integer accumulators; a not-applicable float is `null`.

```
{
"format_version": 2,
"meta": {"samstat_version": "2.4.0", "input_file": "<basename>"},
"summary": {
  "total_reads", "read1_reads", "read2_reads", "paired", "properly_paired",
  "secondary", "supplementary", "is_aligned", "is_partial_report",
  "total_bases", "min_length", "max_length", "n50", "median_length",
  "mean_q", "pct_mapped", "pct_mapq_ge30", "pct_clipped", "median_identity",
  "layout",                      "short" | "long"
  "reads_without_md", "identity_reads", "identity_na", "bases_beyond_position_axis",
  "mapq_counts": {"read1": {<4 ids>}, "read2": {<4 ids>}}
},
"mapq_bins": [{"id": "mapq_ge30", "label": "MAPQ >= 30"}, ...],   ids frozen
"axes": {
  "length_edges": [int],         lower edge of each storage length bin (1922)
  "absolute": {"bins": 1000, "label": "position in read (bp)"},
  "relative": {"bins": 100,  "label": "position in read (% of length)"},
  "identity_bins": 1001,
  "quality_error_scale": 1073741824
},
"read1": {
  "n_reads": int,
  "length":      {"reads": {<id>: [int]}, "bases": {<id>: [int]}},
  "read_mean_q": {<id>: [int]},                      index = Q 0..93
  "identity":    {<id>: [int]},                      index = 0.1% bin; [] for unmapped
  "absolute":    <position block>,
  "relative":    <position block>
},
"read2": { same as read1 }                           present only if read2_reads > 0
}

<position block> = {
  "composition":    {<id>: {"A": [int], "C": [int], "G": [int], "T": [int], "N": [int]}},
  "qual_bases":     {<id>: [int]},
  "qual_error_sum": {<id>: [int]},
  "aligned":        {<id>: [int]},
  "aligned_md":     {<id>: [int]},
  "mismatches":     {<id>: [int]},
  "mismatch_base":  {<id>: {"A": [int], "C": [int], "G": [int], "T": [int], "N": [int]}},
  "insertions":     {<id>: [int]},
  "deletions":      {<id>: [int]}
}
```

- Added after 2.4.0 (optional field, `format_version` stays 2 per `docs/json-format.md`): `mismatch_base`, the per-read-base split of `mismatches`, written after `mismatches`; same trimming rules; the unmapped class has empty arrays.
- Removed from v1: `read*.base_quality` (Q x position matrix), clip-trimmed `length`. Changelog gets one line for v1 consumers.
- Array length bounds (tested by `longread:ultralong`): absolute-axis arrays <= 1000, relative-axis arrays <= 100, identity arrays <= 1001, `read_mean_q` <= 94.
- Size: a 1 Mb single read gives well under 1 MB (about 15k numbers). A 30x WGS file with all 8 class/mate groups full is about 1 MB.
- Transitional (waves 3-4 only, `work-packages.md`): the v2 blocks are first emitted under a top-level `"v2"` key next to the v1 keys, then promoted in WP4.3.

## 8. MultiQC custom content (`--multiqc`)

File names unchanged: `<input>.samstat_<section>_mqc.json`.

**General statistics (7 columns; the 4 existing ids are never renamed; `samstat_n50` is the display value of 3.1 and is left out when no read has a sequence):** `samstat_total_reads`, `samstat_pct_mapq_ge30`, `samstat_pct_unmapped`, `samstat_pct_proper_pair` (paired only), `samstat_n50`, `samstat_mean_q` (omitted for FASTA), `samstat_median_identity` (aligned only).

**Sections (5):**

| id | type | content |
|---|---|---|
| `samstat_general` | generalstats | the columns above |
| `samstat_mapq` | table | reads per MAPQ class, R1+R2 |
| `samstat_length` | linegraph | R1 reads per display length bin (3.1), x = bin lower edge, log x |
| `samstat_basequal` | linegraph | R1 mean Q per bin on the **relative axis** (100 bins), all classes combined |
| `samstat_error` | linegraph | R1 (mismatch + insertion + deletion) per aligned base, %, MAPQ >= 30 only, **relative axis**; aligned input only |

The two position linegraphs always use the relative axis, labelled `position in read (% of read length)`, for every sample whatever its HTML layout, so short and long samples in one cohort are comparable. No per-MAPQ plots, composition or identity histograms in MultiQC; those stay in the per-sample HTML.

## 9. What is not in the report (by design)

GC content, duplication, overrepresented sequences, adapter content, k-mers, per-tile, 2D length x quality or length x identity, coverage, insert size, homopolymers, modified bases, configurable MAPQ classes, multi-sample views (MultiQC does that), pass/warn/fail flags, `--plotend` plots, the MAPQ stacked bar, per-plot data tables, hover tooltips, a first-100-bp inset for long reads (`decisions.md` A7).
