# Topic design: long-read support (SAMStat 2.4)

Status: FINAL design for implementation. Base: the "minimal" proposal. Grafted from the "researcher" proposal: aligned-base denominators so error profiles can be shown as rates (GitHub issue #12), 1% length bins, the rule printed in the report, `n_identity_reads`, a raised Q cap (as a decision), and the merge test.

The guiding rule is that every accumulator has a fixed size, whatever the read length. There are three new fixed-size histograms, a few summary numbers and one automatic switch: the x-axis of the position plots. There is no new option, no mode flag and no new plot family apart from per-read identity.

## 0. Verified code facts this design relies on

All of these were checked against the `planning/v2.4` tree on 2026-10-14.

- **Matrix growth.** `collect_stats` resizes every position matrix to `sb->max_len+1`. This covers `base_comp`, `qual_comp`, `mis`, `ins`, `del`, `len_dist` and the `_end` matrices (`src/collect/collect.c:92-116`). This growth is the source of the memory problem and the ~100 MB JSON.
- **Length histogram.** It stores the clip-trimmed length, `itm->len - (clip_start + clip_end)` (`collect.c:403`).
- **Secondary and supplementary records.** They never reach the collector. `read_bam_chunk` counts them and skips them (`src/htsinterface/htsglue.c:82-86`). That means every histogram is already primary-only, and `n_secondary` and `n_supplementary` already exist (`collect.h:67-68`, `json_report.c:93-96`).
- **`a->error` (NM) is stale or uninitialised.** It is assigned only when an NM tag is present (`htsglue.c:165-169`). Buffer slots are reused, so a record without NM inherits the previous record's value.
- **Mismatches are lost silently without MD.** `parse_alignment` sets `genome[rp] = read[rp]` for M, `=` and X (`src/sambamparse/sam_bam_parse.c:158-166`). Only the MD walk adds mismatches. A BAM without MD therefore shows a mismatch profile of zero, even when the CIGAR has X operations.
- **Insertions and deletions are counted as events.** One count per gap opening (`collect.c:605-621`). This already matches gap-compressed identity.
- **`--norm` exists only for composition and quality.** The 100-bin relative axis is allocated only when `--norm` is set (`collect.c:689-750`, `collect_stats_norm_seqqual` at `collect.c:441`). The quality matrices have 42 rows (`QUAL_MAX_BIN 41`, `collect.c:43`).
- **`config_stats` (`collect.c:248`) is the existing one-shot, first-chunk detection hook.**
- **Chunk size.** `buffer_size` defaults to 10,000 records (`src/param/param.c:390`).

## 1. Decisions

### 1.1 Position axis: absolute for short reads, relative for long reads, chosen once per file

- In `config_stats`, compute the median `l_qseq` of the first chunk. If it is at least 1000 bp, the file is `long`, otherwise `short`. The choice is frozen for the whole file.
- **Short mode.**
  - Absolute positions, with the matrices fixed at 1000 columns.
  - Bases past position 1000 are not counted in the position plots, and the caption says "first 1000 bases".
  - Render-time `clu_data` binning is unchanged, so short-read position output is byte-identical for reads up to 1000 bp.
- **Long mode.**
  - Only the 100-bin relative axis is collected. This is the existing `--norm` path, extended to mismatch, insertion, deletion and the new aligned-base row.
  - The plot titles are unchanged. The x-label is "position along read (% of aligned length)".
  - The bin index is computed per bin boundary (`ceil(k*len/100)`), not per base with a division.
- Only the axis in use is collected, so short-read speed is unaffected.
- The report header prints the rule as one line, for example `Layout: relative positions (median read length 12.4 kb >= 1 kb)`. The JSON carries `"layout": "short" | "long"`.
- **Memory.** It is constant: at most about 3.5 MB per thread (about 440 rows × 1000 columns × 8 B for R1 and R2). A 1 Mb read no longer costs about 3.5 GB.
- **Why not the researcher's geometric `pos_bin`?** It needs a new bin function, new axis labels and new JSON bin edges. The two axes above already exist and cover both regimes.
- **Interface to the I/O and threading topic.** The mode must be decided from the first chunk before work is dispatched to workers. Every worker must use the same mode.

### 1.2 Error profiles as rates (grafted; answers issue #12)

- Add one row per MAPQ bin to the mismatch matrix: aligned read bases at each position (M, `=`, X, I). It is incremented in the existing walk in `get_aln_errors`, which adds no extra pass.
- The report can then show mismatch, insertion and deletion events per aligned base. That rate is comparable between MAPQ bins and between short and long reads, where raw counts are not.
- Whether the three profiles are merged into one plot is decided by the presentation topic. This topic only supplies the denominator.

### 1.3 Read-length histogram: log bins, bases per bin, full length

- **Bins.** Exact (1 bp) below 1000 bp, then 1% geometric bins up to 10 Mb. Lengths of 10 Mb or more go into a final overflow bin. That is about 1926 bins.
- **Contents.** Each bin stores a read count and a base sum, per MAPQ bin, for R1 and R2. That is about 250 KB in total.
- **Length definition.** Full `l_qseq` (`itm->len`), not the clip-trimmed length. N50 is a yield statistic. Records are already primary-only (§0).
- **Derived values.**
  - N50 and mean length come from the bin sums. They are exact below 1 kb and within 1% above. The report shows `~` only when N50 is above 1 kb.
  - Exact scalars are kept next to the histogram: `n_reads`, `n_bases`, `min_len` and `max_len`.
- **Plot.** The x-axis is always log, with N50 drawn as a vertical line. There is no switch: a 150 bp library shows as one honest spike.
- **Short reads.** Unclipped short reads keep exact lengths. Golden values change only for soft-clipped reads.

### 1.4 Per-read identity (the one new plot): CIGAR plus NM, no MD needed

- **Definition.** Gap-compressed identity, `1 - (X + I_ev + D_ev) / (M + I_ev + D_ev)`.
  - M is the number of M, `=` and X bases.
  - `X = NM - I_bases - D_bases`.
  - `N` (intron) operations are ignored.
  - This is minimap2's `de` convention (https://lh3.github.io/minimap2/minimap2.html).
- **Where it is computed.** In the first CIGAR loop of `parse_alignment` (`sam_bam_parse.c:72-100`). The cost is O(CIGAR operations), with no per-base work.
- **Source, in order.**
  1. NM.
  2. If there is no NM but the CIGAR has `=`/X operations, use the X bases directly. This costs nothing.
  3. Otherwise the record is counted in `n_identity_na`.
  - If `X` comes out negative, the NM is inconsistent with the CIGAR. The record is counted in `n_identity_na` and not binned.
  - The `de` tag is not read and there is no MD fallback. That gives one definition and one code path.
- **Histogram.** 1001 bins in 0.1% steps, per MAPQ bin, for R1 and R2, primary mapped records only. That is about 64 KB.
- **Headline.** Median identity comes from the histogram (±0.05%).
- **Missing tags.**
  - When `n_identity_na > 0`, the caption reads "identity unavailable for N of M mapped reads (no NM tag)".
  - If no record has a usable tag, the plot is replaced by that sentence. There is never a silent zero.
- **Fix.** Reset `a->error = -1` per record in `read_bam_chunk` before the NM lookup.
- **Short reads.** The plot appears for short-read aligned input too. One layout serves both.

### 1.5 Per-read mean quality: error-probability space, integer bins

- **Computation.**
  - Per read, take `-10·log10(mean(10^(-q/10)))` over the bases with quality.
  - It is computed in the existing per-base quality loop (`collect.c:417-436`), with a 256-entry lookup table.
  - Records with missing QUAL (0xFF) are skipped.
- **Histogram.** 94 integer bins (Q0–Q93, the SAM maximum), per MAPQ bin, for R1 and R2.
- **Headline.** Median read Q comes from the histogram.
- **Determinism.** A summary mean summed in floating point is deliberately not used. A float sum depends on merge order and would break byte-identical output across `-t`. All accumulators stay integer.
- **Tags.** The `qs` and `rq` tags are not read. Their semantics are UNVERIFIED: the Dorado README does not document `qs` (https://github.com/nanoporetech/dorado). QUAL works for both FASTQ and BAM.

### 1.6 Missing MD: make it visible, and use `=`/X when present

- **Fix.** In `parse_alignment`, set `genome[rp] = 'N'` for X operations, so `--eqx` CIGARs produce mismatches without MD. MD still overrides when present.
- **Count.** Count mapped records that have neither MD nor X operations as `n_no_md`. The error-profile caption shows "mismatch profile excludes N reads without MD tag".
- **Collection.** Mismatch counts are not collected for those reads. Their aligned bases are also left out of the denominator, so the rates stay correct.
- **Why.** Today this case is silent miscounting, which the user-needs research ranks as a top annoyance. CRAM decode of MD (`CRAM_OPT_DECODE_MD`) belongs to the I/O topic.

### 1.7 Supplementary and secondary: two numbers, no plot

- Report `supplementary records per 100 primary mapped reads` and `secondary records per 100 primary mapped reads` in the summary strip and the JSON. The counts already exist.
- There is no SA-tag parsing and no per-read chimera call.

### 1.8 Chunk buffer bounded by bases

- `read_bam_chunk` stops at 10,000 records or 32 Mbases, whichever comes first. Today, 10,000 slots at 100 kb each hold about 4 GB.
- The FASTA/FASTQ reader gets the same budget, or the I/O topic's htslib reader inherits it.

### 1.9 Summary strip and bounded output

- **Summary strip** (the same cells every run; cells that do not apply show "–"):
  - reads
  - bases
  - N50
  - median length
  - median read Q
  - mapped %
  - median identity
  - supplementary per 100
- **New JSON fields.** All are fixed-size arrays with explicit bin edges, so consumers do not re-implement bin functions.
  - `layout`, `n50`, `length{edges,reads,bases}`, `identity{…, definition:"gap-compressed (CIGAR+NM)", n_identity_reads, n_identity_na}`, `read_mean_q`, `n_no_md`, `supplementary_per_100`, `secondary_per_100`.
  - The `schema_version` minor version is bumped. Zero tails are trimmed as now.
- **Size bound.** A 1 Mb read gives a JSON under 1 MB (today about 96 MB).
- **MultiQC.** `_mqc` gains stable general-stats columns: N50, median identity and median read Q.

## 2. Explicitly OUT

- `--long` or any new option. `--norm` becomes automatic. The CLI topic should make it a deprecated no-op for one release.
- Geometric position bins, and absolute positions in long mode.
- Reading the `de`, `qs`, `rq` and SA tags. An MD-based identity fallback. Three identity methods.
- Scatter, 2D and kde plots (length × quality, length × identity), aligned-fraction and soft-clip plots, homopolymer analysis, yield over time, channel maps.
- Pass/fail flags for long-read metrics. Typical-range captions belong to the presentation topic and stay UNVERIFIED until checked on real Dorado and Revio BAMs.
- Lazy R2 allocation (memory is already bounded). Parse-loop speed rewrites (the performance topic).
- Modified bases (separate topic).

## 3. Contract for other topics

- **Threading.** `stat_collection` stays purely integer and additive. Merging is element-wise addition, so output is identical for any `-t`. The mode is decided once, before dispatch.
- **Presentation.**
  - It owns plot rendering, rate display, captions and the render-binning change (sum vs rounded mean).
  - It must not need any data beyond §1.
  - The `_end` matrices also grow with `max_len`. They must either be capped at 1000 columns or removed with `--plotend` (see DECISION NEEDED 4).
- **I/O.** It owns the base-bounded chunk in the htslib FASTQ reader and `CRAM_OPT_DECODE_MD`.

## 4. DECISION NEEDED

1. **Q cap for the per-position quality matrix.**
   - Keep Q41: no golden churn, but HiFi and AVITI saturate the per-position plot.
   - Raise to Q60 with an overflow bin: about 1.5x the rows, and goldens change only for data above Q41.
   - Recommendation: raise to 60. The per-read Q histogram covers Q0–93 either way.
2. **Layout threshold.** 1000 bp, or 2000 bp, on the first chunk's median. Libraries from 500 bp to 2 kb (amplicons, some HiFi) sit near the boundary. Recommendation: 1000.
3. **Start-of-read artefacts in long mode.** At 1% bins, the first bin of a 10 kb read is 100 bp, so adapter signal is diluted. Accept this, or also collect the first 100 absolute positions in long mode (one more axis and plot variant). Recommendation: accept, for simplicity.
4. **`--plotend` and the `_end` matrices.** Remove them, since the relative axis shows read ends, or keep them capped at 1000 columns. This is shared with the CLI topic. Recommendation: remove.
5. **Length change for clipped reads.** Using full length instead of the clip-trimmed length changes 2.3 golden values and the meaning of `len_dist`. Confirm.

## 5. Affected files

- `src/collect/collect.{c,h}`: fixed widths (remove `plot_data_resize_len` growth), mode detection in `config_stats`, the relative axis for mis/ins/del plus the aligned-base row, and the length, identity and read-Q histograms.
- `src/sambamparse/sam_bam_parse.c`: identity counters in the CIGAR loop, X operations as mismatches, `n_no_md`.
- `src/htsinterface/htsglue.c`: reset `a->error`, chunk budget by bases.
- `src/plot/*`: fixed widths.
- `src/report/{stat_report,json_report,mqc_report}.c`: the summary strip, the layout line, and the new plots and fields.
- `docs/samstat.schema.json`, `docs/json-format.md`, `docs/multiqc.md`.

## 6. Effort

About 7 days for one developer:

| Part | Days |
|---|---|
| Fixed matrices, mode switch, relative error profiles | 2.5 |
| Identity, X/MD handling, NM reset | 0.75 |
| Length and read-Q histograms | 0.75 |
| Report, JSON, MQC and schema | 1.5 |
| Fixtures, tests, one golden-regeneration commit | 1.5 |

Order: bounded model first (it fixes the out-of-memory failure), then length/N50, then identity.

## 7. Risks

- **Golden churn** from the full-length change and the new fields. Regenerate in one dedicated commit.
- **Mixed-length files.** If the first chunk is short but the file holds many long reads, bases past 1 kb are lost from the position plots. Length, N50 and identity stay correct, and the caption says so.
- **NM availability in pbmm2 and Dorado BAMs** is UNVERIFIED. The n/a path covers it. Check a real Dorado R10.4.1 BAM and a real Revio BAM before release.
- **Identity differs from cramino and NanoPlot**, which use different definitions. The formula is printed and named in the JSON.
- **1% length bins** make N50 approximate above 1 kb (by at most 1%). The exact scalars are kept anyway.

## 8. Tests

1. **Fixtures.** A seeded generator (`uv run python`) writes fixtures to `tests/data/longread/`, under 500 KB in total. Expected values come from the generator's ground truth, never from samstat output.
   - An ONT-like BAM: 1–100 kb reads with known identity, with and without NM.
   - An `--eqx` BAM without MD.
   - A supplementary and secondary chain.
   - A long FASTQ.
   - A single 1 Mb read.
2. **Accuracy.**
   - N50 is exact below 1 kb and within 1% above.
   - Median identity is within 0.1% of the truth.
   - Supplementary and secondary counts are exact.
   - The X-only fixture gives a non-zero mismatch profile.
3. **Bounds.** A 1 Mb read gives JSON and HTML under 1 MB and peak RSS under 200 MB (Linux, `perf` label, ratio-gated).
4. **Layout.** A 150 bp file reports `layout: short` and a 20 kb file reports `long`.
5. **Short-read invariance.** For unclipped short-read goldens, position arrays are byte-identical. The diff is limited to the new fields and to `len_dist` for clipped reads.
6. **Missing tags.**
   - With no NM, the caption shows identity n/a.
   - With no MD and no X, the caption shows the `n_no_md` count.
   - In both cases, nothing crashes under ASAN, and the stale `a->error` is gone.
7. **Merge.** Merging two collections equals serial processing (shared with threading).
8. **Report content.** Every report contains the fixed plot titles. `check_json.cmake` needs `>=` and `<=` operators (owned by the test topic).
