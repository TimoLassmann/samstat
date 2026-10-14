# SAMStat 2.4 plan review log

Date: 2026-10-14. 55 issues from four review lenses (feasibility 17, performance 9, researcher 13, scope 16). Verdicts: **A** accepted, **P** partly accepted, **R** rejected. Files: SP `standard-plots.md`, WP `work-packages.md`, DE `decisions.md`, RM `README.md`. Spot checks against the code: `samstat.c:178` (`is_aligned = 1`), `stat_report.c:52/315` (file stat), `collect.c` `plot_data` use and `get_mapqual_bins` (MAPQ >= 30 includes 255), `run_collect.cmake:38` / `run_json.cmake:54` (`format_version 1`), `htslib.cmake:25` (`BUILD_COMMAND make`): all as the reviewers said. The length grids were computed (`uv run python`): storage 1922 bins; `g(1)` does not pass through 1000; a single 1%-from-100 grid puts 150 bp in bin [149,151).

## Feasibility

| # | Issue | V | Change |
|---|---|---|---|
| F1 | WP4.3 deletes v1 keys but does not own tests that read them | A | WP4.3 owns `tests/{collect,reader,fastx}/CMakeLists.txt`, `run_collect.cmake`, `run_json.cmake`, `check_schema_keys.cmake`; Do item to rewrite or delete every v1 CHECKS path and list them in the PR (WP). |
| F2 | FASTQ via htslib would look aligned, /2 reads to R2 | A | New SP 3.0 input contract (no `aln_data`, flags ignored, all R1); WP3.2 implements it and removes `is_aligned = 1`; `fastx:r1r2_suffix` enabled in WP2.1 to pin 2.3 behaviour. |
| F3 | `cli:stdin_*` fail on `stat(\"-\")` in `stat_report.c` | A | WP3.2 owns one hunk of `stat_report.c` (skip file-stat section for `-`). |
| F4 | `plot_data` used by collect/json/mqc; WP4.2 could break it, WP5.1 cannot delete it | A | WP4.2b only adds files; `plot_data` unchanged until WP5.1, which now owns `src/plot/plot.*`, `plot_group.*`, `plot_interval.*`, `plot_test.c`, `src/plot/CMakeLists.txt`. |
| F5 | WP5.1 needs `samstat.c`/`param` owned by WP5.3 | A (option b) | WP5.1 runs after WP5.3 and owns `samstat.c` and `param.{c,h}`; rule 1 allows sequential "after X" packages. WP3.1a keeps `stat_collection_config_additional_plot` until then. |
| F6 | Footer command line not stored; breaks `-t` HTML identity | A (simplest form) | Footer has no command line (version, basename, citation, counts) (SP 1, DE F10). No `p->cmdline`. |
| F7 | Display and storage length grids share no edges above 1000 | A | One master grid: display = exact 1-100, 1% to 996, then storage edges from 1000; bin count corrected to 1922; unit test in WP3.1a (SP 3.1, DE F2). |
| F8 | Content detection vs empty-input tests | A | Empty/whitespace input = exit 0, total 0, before detection; empty-sequence records not counted (SP 3.0, WP3.2 acceptance lists the four tests). |
| F9 | `is_aligned` undefined | A | Defined in SP 3.6: `has_sq` or any primary mapped; fixtures `ubam_like.sam`, `sq_nomap.sam` (WP1.2), asserted in `report:na_notes` and `longread:*`. |
| F10 | WP3.1 and WP4.2 are XL on the critical path | A | Split into WP3.1a (collector) / WP3.1b (JSON preview, after 3.1a) and WP4.2a (page) / WP4.2b (plot renderers, fixed interface, parallel). |
| F11 | `config_stats` first-chunk logic makes v1 depend on batches | A | Moved into WP3.1a: decisions per record (SP 3.0). |
| F12 | Ultralong "same array lengths" ill-defined; SAM generation unspecified | A | Bounds (<=1000/100/1001) and `string(REPEAT)` in `run_ultralong.cmake` (WP2.1, SP 7). |
| F13 | No XML tool for `svg:xml_wellformed` | A (alternative) | Tag-balance checker inside `svg_plot_test.c`; no xmllint, no CI change; test folded into `svg:unit`. |
| F14 | Mean Q definition: all bases vs relative-axis sums | A | Merged with S5: scalar pair over all bases incl. clips; plots 3-5 aligned part only (SP 1.2, 3.3, DE F12). |
| F15 | `test_view` not built by default for bench | A | `bench.sh` uses SAM input; generation time recorded separately (WP2.3). Small fixture BAMs keep the existing `tests/data/README.md` method. |
| F16 | Protein FASTA is not rejected, only garbled | P | Wording fixed (DE B1: unsupported, shows near-100% N). Rejected the exit-1 detection: htslib has already mapped letters to N, so detection needs extra raw-text parsing for an input nobody runs on this tool. |

## Performance

| # | Issue | V | Change |
|---|---|---|---|
| P1 | Thread budget exceeds `-t`; absolute targets not credible | A | Formula in WP4.1 (T=4: reader, 1 pool, 2 workers; T=8: P 2, W 5; n_files >= T keeps 2.3 scheduling). Targets as ratios on x86-64 Linux in WP2.3; prototype table not quoted (RM goal 3, DE B7). |
| P2 | Merge of non-additive state | A | Per-field merge semantics in SP 3.6; `min_length` = UINT64_MAX; `is_aligned` derived in finalise; zero-read merge unit test (WP3.1a). |
| P3 | Thread test fits in one batch | A | Undocumented test-only env vars `SAMSTAT_BATCH_RECORDS/BASES` (WP3.2); `cli:threads_one_file` uses them and asserts `workers: N` > 1 (WP2.1, WP4.1). TSAN runs this test. |
| P4 | Per-base cost of v2 collector | P | Incremental relative-bin walk, binary search for Q bin, array of structs (SP 3.2, 3.3). Gate is relative and manual: <= 25% slower than 2.3 at `-t 1`, reported in the PR. No ns/base micro-benchmark harness. |
| P5 | FASTQ.gz target; build differences | A | Relative target (>= 1.5x 2.3.0), x86-64 Linux; `bench.sh` prints htslib version and libdeflate; `docs/benchmark.md` records the build. Moving FASTQ parsing out of the reader not adopted (measure first). |
| P6 | Memory goal has no automated gate | P | Queue = W+2 batches and the memory bound written in WP4.1; 8 x 1 Mb case with a 1 GB target in `bench.sh`. CTest RSS check rejected (platform-specific, breaks under ASAN; DE E2 left to the maintainer). |
| P7 | HTML size worst case is short paired input | A | Merged with S1: no per-point tooltips; x range to last non-empty bin; `report:short_html_size` on a paired 999 bp file, < 1 MB. Point thinning rejected (not needed without tooltips: about 36k points x about 12 B). |
| P8 | uint64 overflow of summary sum | A | `qual_error_sum_all` is `unsigned __int128`, internal only; unit test above UINT64_MAX (SP 3.3, WP3.1a). |
| P9 | `-p` counting must stay in the reader | A | Reader-only counting and mid-batch stop (WP3.2, WP4.1); `reader:p_cutoff_counts` at `-t 1` and `-t 8` in WP4.1 acceptance. |

## Researcher

| # | Issue | V | Change |
|---|---|---|---|
| R1 | Plot 5 empty for most BAMs (no MD) | P | Note and N/A text give `samtools calmd` and `minimap2 --MD` commands (SP 1.3, plot 5); note already in the Notes block. `--reference` stays out of 2.4, but G5 counts real BAMs without MD and DE B2 is flagged for revisit. |
| R2 | uBAM looks like a failed mapping | A | Merged with F9; `@SQ` with zero mapped keeps aligned layout and shows `No mapped reads.` for plots 5-6 (no new note). |
| R3 | Log-x length axis bad for short reads | A | Linear x in short layout, log in long; minimum span 10 bp / one decade; `single_len.sam` asserted (SP plot 1). |
| R4 | Identity comb and 100% spike for short reads | P | 1% display bins in short layout (SP plot 6). Separate labelled 100% bar and y clipping rejected (extra rendering rules; caption already says short reads pile up at 100%). |
| R5 | Plot 4 colours reuse MAPQ colours | A | Disjoint palette (A green, C reddish purple, G black, T brown, N sky blue dashed), swatch row, palette-disjoint unit test (SP plot 4, WP4.2b). Direct A/C/G/T/N labels replaced by the swatch row (see S2). |
| R6 | Error-probability mean vs FastQC | P | Captions of plots 2 and 3 say it reads lower than FastQC; FastQC column says "different averaging". `docs/interpreting.md` with example reports rejected: `simplicity-review.md` deferred a second copy of the help text. |
| R7 | MultiQC mixes bp and % axes | A | Merged with S4: relative axis for every sample (SP 8, DE F13). |
| R8 | MAPQ wording and 255 | A | Benign line reworded; MAPQ caption sentence; 255 classed >= 30 as in 2.3 (`get_mapqual_bins`), `mapq255.sam` fixture and test (SP 1.4, 5). |
| R9 | Clipping invisible; median length duplicated | A | Strip slot 4 = `% bases clipped` (aligned only), median length stays in header and JSON; plot 6 caption says clips not counted (SP 1.2, DE F14). |
| R10 | Long layout hides read starts; silent switch | P | Header shows the reason (`long when >= 1 kb`); axis labels already say "% of length"; plot 3 help text names the limit. Inset panel rejected (DE A7). |
| R11 | Direct labels collide | R | Moot: direct labels removed (S2). |
| R12 | Reads vs pairs | P | Label `Reads (R1+R2)` on paired input. Extra "pairs" sub-value rejected (one number per slot). |
| R13 | Read-only input directory | A | Error message names `-d DIR`; `cli:input_dir_readonly` test (WP1.4, SP 6). |

## Scope

| # | Issue | V | Change |
|---|---|---|---|
| S1 | Per-point tooltips are interactivity and bloat | A | Removed; one series-level `<title>` (SP 2.6, WP1.6, DE F11); size test added (P7). |
| S2 | Two legends; label collision code | A | No direct labels; MAPQ table is the key; plot 4 swatch row (SP 2.4, WP1.6). |
| S3 | Unaligned layout contradicts across files | A | Unaligned = plots 1-4 only, no N/A lines; N/A only for applicable plots lacking data (SP 1, 4; WP4.2a `report:na_notes`). |
| S4 | MultiQC mixed axes | A | Relative axis always (SP 8; WP4.3 test `multiqc:relative_axis`). |
| S5 | Mean Q and clipping definitions disagree; adapter hint wrong | A | See F14; plot 4 hint limited to unaligned input. |
| S6 | Short-layout x range unspecified | A | x = 1 to last bin with any base, shared R1/R2; asserted on `single_len.sam` (SP 3.2). |
| S7 | Plot 3 averaging not explained | A | See R6. |
| S8 | CRAM is not free | A | Cost stated in WP5.2; first package dropped on slip; DE B8 for the maintainer. Not deferred outright: CRAM is a common cohort format. |
| S9 | `-t` total threads false; W undefined | A | See P1. |
| S10 | README line 43 contradicts A3 | A | Reworded (RM). |
| S11 | Plot 5 R1/R2 y-scale | A | Each row its own scale, shared between R1 and R2 (SP plot 5). |
| S12 | Use one length grid | R | A 1%-from-100 grid reports a 150 bp library as N50 149 (computed). F7's master grid keeps exact lengths below 1 kb with no extra frozen surface beyond `length_edges`. |
| S13 | Note texts differ from simplicity-review | A | SP 1.3 declared the single source; precedence line in RM and SP header. `simplicity-review.md` itself not edited (out of this task's file list). |
| S14 | Modbases non-goal vs roadmap | A | Non-goal "in 2.4"; removed from the 2.5 row and WP §5 (RM, WP). |
| S15 | Workflow sized to agent count | A | Opus reviewers only on M/L or Opus packages; 42 runs (WP §1, DE F9). |
| S16 | TSAN job flakiness | A | Limited to `-L golden` and `-R threads`, sysctl step (UNVERIFIED issue link) (WP2.3). |

## Not changed

`simplicity-review.md`, `topics/` and `research/` were not edited. Where they now disagree (simplicity-review §2 strip list with median length, §3 C1 and C18 MultiQC axis, §4 note texts, §5 CRAM "free"), the precedence line in `README.md` and the entries above decide.
