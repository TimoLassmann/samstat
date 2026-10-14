# code-collect: statistics data model and parsing in SAMStat 2.3.0

Scope: `src/collect/collect.{c,h}`, `src/sambamparse/sam_bam_parse.{c,h}`, plus the glue they depend on (`src/htsinterface/htsglue.{c,h}`, `src/plot/plot.c` `plot_data`, `src/plot/plot_interval.c`, driver loop `src/samstat.c`). Paths relative to the repo root; line numbers are for the tree at branch `planning/v2.4`. No timings were measured here; every performance statement is from reading the code and is marked as an estimate where numeric.

## 1. Data model

### 1.1 One matrix per plot, position axis = absolute read position

`struct stat_collection` (`src/collect/collect.h:23-75`) is a bag of `struct plot_data*`, one per plot and per mate (R1/R2): base_comp, qual_comp, mis, ins, del, len_dist (`collect.h:25-41`), optional `*_end` (`:44-48`) and `*_len_norm` (`:51-55`). Plus `basic_nums[2][4]` (counts per mate x MAPQ bin, `:61`) and a few `uint64_t` counters (`:62-68`).

`struct plot_data` (`src/plot/plot.h:15-43`) is both the accumulator and the render config: `uint64_t** data` of `L` rows x `len` columns, plus label strbufs, plus render-time `clu_data`. So the collector and the plotter share one struct; the data model cannot be changed without touching `plot.c`.

Rows (`L`) per plot, from `stat_collection_alloc` (`collect.c:906-928`), with 4 MAPQ bins (`get_mapqual_bins`, `collect.c:1234-1270`: >=30, 10-29, <10, unmapped; `n_bin = 4` at `:1243`):

| plot | rows | meaning | alloc line |
|---|---|---|---|
| base_comp | 5 x 4 = 20 | A,C,G,T,N per MAPQ bin | `:907` |
| qual_comp | 42 x 4 = 168 | Q0..Q41 per MAPQ bin | `:911` |
| mis | 5 x 4 = 20 | mismatching *read* base per bin | `:915` |
| ins, del | 4 each | one row per MAPQ bin | `:919,923` |
| len_dist | 4 | per MAPQ bin | `:927` |

Total 220 rows per mate, 440 rows for R1+R2, each cell a `uint64_t`: **3520 bytes per read position** for the standard plots. Row addressing is `data + group_size * q_idx` then `[base-or-Q][pos]` (`collect.c:301,331,410,423`).

Optional plots (`--end`, `--norm`) are allocated on top: end plots fixed at 250 columns (`:762-767`), normalised plots fixed at 100 columns (`:694-698`). These are already bounded and length-independent.

### 1.2 Growth with read length: unbounded

`collect_stats` resizes every standard plot to `sb->max_len + 1` columns on every chunk (`collect.c:92-108`; `plot_data_resize_len`, `src/plot/plot.c:1072-1095`, allocates a new matrix, copies, frees; it does nothing if already large enough, `:1074`). `sb->max_len` is raised by `read_bam_chunk` (`htsglue.c:191-193`) and is never reset (`reset_tl_seq_buffer`, `src/tld/src/seq/tld-seqbuffer.c:317-332`, does not touch it). Consequences:

- Memory is `3520 B x longest read seen`. 150 bp: ~0.9 MB. 10 kb: 35 MB. 100 kb: 352 MB. **1 Mb: ~3.5 GB per thread**, allocated the moment one such read appears, even in a file where it is the only long read. `-t` threads multiply this (one `stat_collection` per file/thread, `samstat.c:170`).
- R2 matrices are resized and allocated even for single-end data (`collect.c:93,96,99,...` unconditional).
- Resize is exact (`max_len+1`), not geometric, so a file whose reads get progressively longer reallocates and copies all 12 matrices at each new maximum (O(max_len) per event; bounded by the number of distinct new maxima, in practice small, but a 1.7 GB copy when it hits).
- Reads beyond the allocated width are silently truncated for composition/quality/mismatch (guards at `collect.c:302,333,417,602`), but since width = max_len this never triggers for standard plots.

### 1.3 Length histogram: one column per base pair

`len_dist_*` is indexed by exact length, `clamp_idx(len - clips, len_dist_R1->len)` (`collect.c:295,326,403`). So the length histogram has `max_len+1` columns x 4 bins. Exact for 150 bp, useless as a plot for nanopore (thousands of 1-bp bins, mostly empty), and the root of the 100 MB JSON seen with a 1 Mb read (`json_report.c:177-190` emits `pd->data` trimmed only after the last non-zero, `json_width` at `:296`). Two further problems:

- The recorded length is `itm->len - (clip_start + clip_end)` (`collect.c:403`), i.e. the **soft-clip-trimmed** length for aligned reads, not the read length. For N50 / "read length distribution" the full read length (`b->core.l_qseq`, `htsglue.c:105`) is what users expect; clipped length is an alignment metric. The two diverge on noisy long-read data.
- `len_dist` has the PLOT_MOD_DENSITY flag (`collect.c:1104`) and bypasses `clu_data` binning (`plot.c:58`), so it is never rebinned at render time.

### 1.4 Render-time binning already exists, but after the memory is spent

`plot_add` calls `clu_data` (`plot.c:58-60`) which bins columns for display: `len <= 250` exact; `<= 1000` linear intervals; otherwise exponential (`plot_interval.c:34-60`, FastQC-style, `bin_start = 250`, `target_n_clu = 500`, `plot.c:1040,1042`). `summarise_data` (`plot_interval.c:67-107`) then replaces each bin by the **rounded mean of the counts** in it (`:96-100`) rather than the sum. That is fine for fractions but is lossy for counts and needs the full `max_len` matrix first. This is the key observation for the redesign: **the binning function exists and is already what users see; it just runs at report time instead of at collection time.** Moving the position->bin map into `collect` (an O(1) table or arithmetic function) bounds memory without changing what the plots look like for short reads (exact below 250 bp today).

## 2. Per-record flow and cost

Driver (`samstat.c:181-217`): `read_bam_chunk` fills a `tl_seq_buffer` of `buffer_size = 10000` records (`param.c:390`), then for each record `parse_alignment` + `fix_orientation` (`:200-203`), then `collect_stats`, then `clear_aln_data` + `reset_tl_seq_buffer`. No `hts_set_threads` / thread pool anywhere in `src` (grep: none), so BGZF decompression is single-threaded as well as everything else per file.

Per record, in order:

1. `htsglue.c:97-176`: copy flag/MAPQ; `tld_append` of the name; **decode 4-bit SEQ to chars in a scratch buffer then `tld_append` it** (`:117-122`, two passes + memcpy); same for QUAL (+33 conversion, `:151-155`), or a `0xFF` fill when QUAL is `*` (`:133-138`); `bam_aux_get` for `NM` and `MD` (each a linear scan of the aux block, `:165-175`); copy the CIGAR (`:185-188`). Scratch buffer grows 1 KB at a time while hunting for the longest read (`:109-116`), O(n^2/1024) reallocs worst case but only on new maxima.
2. `parse_alignment` (`sam_bam_parse.c:50-256`): pass 1 over the CIGAR for `q_len` and `aln_len` (`:72-100`), sanity check `q_len == s->len` (`:102`; failing records are treated as unaligned without warning), soft-clip detection (`:107-123`), then **materialise two byte arrays `read[]` and `genome[]` of alignment length** (`:143-186`, 2 bytes per aligned column, `galloc` re-allocation to 1.5x growth `:127-135`), then walk the MD string patching `genome[]` (`:196-231`, `isdigit` per char; MD absent -> `genome = read`, so mismatches silently count as zero, `:191-195`).
3. `fix_orientation` (`:258-277`): reverse-complement `seq`, reverse `qual`, reverse-complement `read[]` and `genome[]` of reverse-strand records, so everything is in read-5'->3' orientation. That is four more full passes for half the records (`tld-revcomp.c:53-62`).
4. `collect_stats_all` (`collect.c:346-439`): per base `t[nuc_to_internal[seq[i]]][idx]++` and `get_base_qual` + `t_q[q][idx]++` (`:417-436`); `get_aln_errors` (`:570-632`) walks `read[]`/`genome[]` again, char compares and `mis/ins/del` increments. Each base touches a different row of a separately allocated row array. For qual, 42 rows x 8 B x `len` => the increment target jumps by `8 x len` bytes depending on Q value; at 150 bp rows fit in L1/L2, at 100 kb+ almost every increment is a cache miss (estimate, UNVERIFIED by measurement).

Rough pass count per aligned base: SEQ decode, SEQ copy, QUAL conv, QUAL copy, CIGAR expand (2 arrays), MD walk, optional reverse (4 arrays), base-count loop, qual loop, error loop: about 10-14 touches of each base, with several heap writes of read-length size. Nothing is SIMD; the work is branchy per base. The per-record cost is linear in length (as it must be), so the redesign target is a constant-factor cut (fewer passes, no `read[]`/`genome[]` expansion) and parallelism, not a better complexity class.

### Memory held by the chunk buffer

Each of the 10000 slots owns growable `seq`, `qual`, `name` strbufs (tld) plus an `aln_data` with `read[]`, `genome[]` (`htsglue.c:278-306`). Slots keep their high-water capacity across chunks (`reset_tl_seq_buffer` only zeroes lengths, `tld-seqbuffer.c:322-326`). With reads of mean length L, every slot trends to hold ~ 4L bytes (seq, qual, read, genome) at its own maximum, so ~ `10000 x 4 x L`: 6 MB at 150 bp, ~400 MB at 10 kb, ~4 GB at 100 kb, with the high-water effect making it worse than the mean suggests. `buffer_size` should depend on read length (for example cap total bases per chunk, ~64 MB) rather than being a constant 10000. (`stat_collection` memory in 1.2 is on top of this.)

## 3. Things in the model that are wrong or missing for long reads

- **No N50 / total bases / min / max / mean**: only the per-bp histogram exists; n_total etc. are computed at `collect.c:658-671` from `basic_nums`. N50 is derivable from the exact histogram but not from anything bounded.
- **No identity or per-read error rate.** `NM` is fetched (`htsglue.c:164-169`) into `a->error`, but nothing reads it (grep of `->error` outside tld: only the write). `alloc_aln_data` (`:278-306`) does not initialise `error`, so for records without an NM tag it is an uninitialised (stale) value. A usable per-read identity is available for free from `NM` and the CIGAR (BLAST-style `1 - NM / (M+I+D)`, gap-compressed variants possible), no MD needed; MD-derived mismatches are only needed for the per-substitution-type profile.
- **No per-read mean quality.** Quality is only a position x Q matrix. Long-read QC conventionally shows a per-read mean-Q histogram (mean taken in error-probability space, 256-entry lookup). Not in the data model.
- **Quality cap Q41** (`QUAL_MAX_BIN 41`, `collect.c:43`): modern ONT/PacBio HiFi go above Q41 (HiFi to Q93 in BAM). Everything above is lumped into Q41 (`:63-65`). Same for Illumina-era assumptions: `base_quality_offset = 33` is hard-coded for BAM (`samstat.c:183-184`), fine for SAM/BAM by spec.
- **Error profile x position for long reads is not meaningful** (as in the task brief): position along a 50 kb read has no biology. The `--norm` plot (`collect.c:441-508`, 100 relative-position bins, `floor(p / t_len * 100)` per base, a `double` multiply and `floor` per base, so roughly an order slower than the plain loop) shows what is wanted for long reads, but only for base composition and quality, not for mismatch/ins/del. Option for the redesign: error rates along the relative position axis for long reads.
- **Errors by position count events, not rates.** `mis/ins/del` columns are raw counts (`collect.c:608-626`); normalising by coverage at that position is left to the plotter. Deletions are attributed to the *next* read position since `rp` is not advanced on `D` (`:614-621`); consecutive deleted columns count once per run (`r[j-1] != '-'` check), i.e. the deletion plot counts deletion *events*, mismatches count *bases*; identity needs bases, so keep both notions explicit in the new model.
- **MAPQ stratification with long reads**: 4 bins x 42 quality rows x positions is the memory driver (168 of 220 rows). For the redesign quality could be stored 4 bins x ~42 rows x ~300 position bins = 400 KB, so stratification stays affordable once positions are binned.
- **Mate assignment**: a record with neither READ1 nor READ2 flag defaults to R1 (`collect.c:355,379-384`), correct for single-end. R2 allocations are not conditional (see 1.2).
- **Dead debug path**: `collect.c:368-376` has an unreachable "Oh no I should reverse!" `LOG_MSG` (`fix_orientation` already clears `reverse`, `sam_bam_parse.c:270`); if it ever fired it would log once per record. Safe to delete during the rewrite.
- **Merge across threads/files**: none exists; each file gets its own `stat_collection` and the report is per file. A bounded model makes a `stat_collection_merge(dst, src)` a handful of element-wise adds (sub-MB), which is what enables within-file parallelism (reader thread + N worker threads each with a private collection, merged at the end).

## 4. Recommended data model (one model for 50 bp and 1 Mb)

Principle: every array size is a compile-time or small constant, independent of the longest read. Everything is a plain sum of `uint64_t`, so merging is addition.

1. **Position axis: a shared bin function `pos_bin(p)`**, exact for `p < 250` (identical to today's output for short reads), then geometric growth (~5% width) up to a hard cap (clamp the last bin at e.g. 10 Mb). That is ~250 + ln(10^7/250)/ln(1.05) ~ 470 bins. Implement as a precomputed `uint16_t` table for p < 64k plus arithmetic (`clz`/`log2` fixed-point) above, or simply table the first 64 K and put everything above in log bins. Store **sums** per bin plus a per-bin positions-per-bin width (known statically) for rate normalisation; drop the rounded-mean step in `summarise_data`. Matrix size: 440 rows x 470 x 8 B = 1.6 MB per thread, constant. Replaces the `plot_data_resize_len` calls (`collect.c:92-116`) and `max_len` dependence; `clu_data` then is only a label generator.
2. **Length histogram: log-binned**, exact below 100 bp then 1% geometric bins (`bin = 100 + floor(ln(len/100) / ln(1.01))`, ~1500 bins to 10 Mb, 12 KB per MAPQ bin), counts **and bases-per-bin** (sum of lengths). From the bases-per-bin array, N50/N90, mean, total bases and a "bases vs length" weighted histogram follow with <=1% error and zero extra passes. Also keep exact scalars `n_reads, n_bases, min_len, max_len`. Use full read length (`l_qseq`), not clip-trimmed; record clipped fraction separately if wanted.
3. **Identity**: per aligned primary record, `ident = 1 - NM/(M+I+D)` from `NM` (fallback: count from the already built `genome`/`read` or from CIGAR+MD), into a fixed 0..100% histogram (1000 bins, 0.1% resolution) x MAPQ bin (x R1/R2). 8 KB per bin. Distribution plot + median/mean in the summary.
4. **Per-read mean quality**: one histogram, Q0..Q60+ in 0.5 steps, mean in error-probability space via a 256-entry table. Raise `QUAL_MAX_BIN` (currently 41) to ~60 or add one "above" row for modern chemistry; this grows qual rows from 42 to ~61 (positions are bounded now, so affordable).
5. **Keep** the MAPQ stratification (distinctive strength) but derive R1/R2 allocation lazily (`n_read2 > 0` on first mate-2 record).
6. **For long reads, error profiles use the relative-position axis** (already in the code as `--norm`, 100 bins) extended to mismatch/ins/del, and the absolute axis only for the first N bins (adapters/ends). Decide once from the first chunk (the code already does one-shot detection in `config_stats`, `collect.c:248-277`): median length > ~1 kb => relative axis default. Alternative with no mode switch: always collect both absolute (bounded, item 1) and relative (100 bins) and let the report pick by length distribution; costs one extra array and one extra `floor`-free index (precompute `idx = (p*100)/len` with integer math, replacing the `double` division per base at `collect.c:494`).
7. **Cost**: (a) drop `read[]`/`genome[]` expansion, and walk CIGAR+MD (or only CIGAR+NM+base compare) directly while counting, one pass over bases; (b) skip `fix_orientation` on the byte arrays by indexing from the other end (`idx = len-1-i`) instead of reversing four buffers; (c) decode SEQ straight from the 4-bit BAM representation into the counters, avoiding the char copy into `tl_seq` for BAM entirely; (d) bound the chunk by total bases, not 10000 records; (e) one `hts_set_threads` / `htsThreadPool` for BGZF decode (standard htslib API; not currently used - UNVERIFIED speedup on the target data, though commonly 2-3x for BAM), then per-chunk worker threads with private collections merged by addition. None of this changes the output for existing short-read goldens except where numbers legitimately change (clip-trimmed vs full length, mean->sum binning).
8. **Row layout**: store the accumulators as one contiguous array with the position index as the **fastest-varying** index of a *small* block (`[bin][pos]` where pos bins <= 500) so a row is <= 4 KB and the working set (1.6 MB) stays in L2; with bounded positions the cache-miss problem in section 2 disappears without further tuning.

## 5. Impact on other modules (for the other planners)

- `plot_data` doubles as the render config; if collect owns fixed-size typed arrays, `stat_report.c`, `json_report.c` (`:177-250` read `pd->data`/`pd->len`) and `mqc_report.c` need to read the new structs; JSON schema (`docs/samstat.schema.json`) changes shape for `length` and position arrays (bin edges become explicit; JSON size becomes bounded: ~0.5 MB regardless of read length).
- Goldens in `tests/golden` will change wherever binning/mean semantic changes; exact-below-250 positions can stay byte-identical for the common short-read tests if the first 250 positions keep the same values.
- FASTA/FASTQ path shares `collect_stats` (`samstat.c:383`) via `collect_stats_seq/seqqual` (`collect.c:279-343`), so the same bounded model must cover it; these paths only fill R1 and UNMAP bin, and index base composition by `itm->len` (not clip-trimmed).
- External facts: only htslib API (`hts_set_threads`, BAM aux `NM`/`MD` tags, SAM spec for tag semantics) are referenced; spec: https://samtools.github.io/hts-specs/SAMv1.pdf and https://samtools.github.io/hts-specs/SAMtags.pdf. HiFi Q93 maximum in BAM: UNVERIFIED here (PacBio documentation not fetched).
