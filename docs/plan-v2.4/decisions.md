# SAMStat 2.4: decisions for the maintainer

Every DECISION NEEDED item from the nine topic files and the simplicity review, deduplicated. Each item: options, recommendation, the ruling `work-packages.md` assumes, and its source. Change a ruling here before gate G0 and the affected packages follow it.

The simplicity review's cuts are binding. **No exception to a cut is taken in this plan.** Section F lists the details this master plan added where the topics were silent; none adds a plot, option or output file.

Legend: **[R]** = ruling applied in the plan. Source tags: SP standard-plots, LR long-reads, PF performance, PR presentation, EC ecosystem, MB modbases, CL cli-ux, TS testing, RC release-cleanup, SR simplicity review.

---

## A. Report content

| # | Question | Options | Recommendation [R] | Sources |
|---|---|---|---|---|
| A1 | Status flags in 2.4? | (a) none, only factual notes; (b) OK/CHECK/FAIL table (max 5 rows) plus `[OK]` on stderr | **(a).** Every proposed threshold is UNVERIFIED in its own topic; uncalibrated red lights are FastQC's known failure (`research/fastqc.md`). Revisit in 2.5 after calibration on real data. | SP D1, PR 3, SR D-2 |
| A2 | Keep base composition? | (a) keep as plot 4; (b) cut to 5 plots | **(a).** It was the paper's headline example (N spike at position 40) and is a FastQC module researchers already read. | SP D2, PR 2 |
| A3 | Base composition: split by MAPQ? | (a) one panel over all reads (per mate); (b) one panel per MAPQ class, as in the paper | ~~(a).~~ **Revised 2026-10-18: (b)**, as short rows (`standard-plots.md` encoding rule 9). Ruling (a) left plot 4 unable to show the 2023 paper's own example (N in confidently vs poorly or unmapped reads), and plot 5 could not carry it either once F3 dropped the base. Rows instead of panels keep the page short. | PR 4, SR D-1, C8, maintainer 2026-10-18 |
| A4 | Paired input: which mates? | (a) R1 and R2 side by side in every plot; (b) R1 only for base composition | **(a).** One rule for all plots (SR encoding rule). Supersedes SP D6. | SP D6, SR §2 |
| A5 | Length for plot 1 and N50 | (a) full read length (`l_qseq`); (b) clip-trimmed (2.3, `collect.c:403`) | **(a).** N50 is a yield statistic. Changes 2.3 values for soft-clipped reads (in golden regeneration #1). | SP D3, LR 5, SR D-5 |
| A6 | Short/long layout switch | (a) median length >= 1 kb, no override; (b) 2 kb; (c) add an override option | **(a).** Decided at end of run (works with stdin). Add an override only after a real misdetection report. | SP D4, LR 2, SR D-6, CL |
| A7 | Start-of-read artefacts in long layout | (a) accept dilution in 1% bins; (b) also collect first 100 absolute positions; (c) an extra first-100-bp inset panel | **(a).** Known limit, stated in the plot 3 help text. The absolute axis is in the JSON anyway (both axes are always collected). No inset: one more panel per plot. | LR 3, plan review |
| A8 | Quality cap | (a) keep Q41; (b) Q60 + overflow; (c) Q93 | **Resolved, no decision left.** The per-position Q histogram is replaced by integer error-probability sums (`standard-plots.md` 3.3); the per-read histogram spans Q0-93. | LR 1, SP §3, SR C4 |
| A9 | Rendering | (a) hand-written inline SVG in C; (b) inlined plotly-basic (about 1.2 MB per report, https://github.com/plotly/plotly.js/blob/master/dist/README.md, via topic research, UNVERIFIED here) | **(a).** Only (a) meets HTML < 1 MB for the ultralong test and needs no JavaScript. Fallback (b) if WP1.6 overruns badly. | SP D5, PR 1, RC 3, SR D-4 |
| A10 | `--plotend` | (a) hidden no-op until 3.0, matrices removed; (b) cap at 1000 columns; (c) make default for short reads | **(a).** The matrices grow with read length; the relative axis shows read ends. | LR 4, CL 1 |
| A11 | Modified bases | (a) defer entirely; (b) two summary rows (reads with MM/ML, reads with inconsistent MM/ML) | **(a).** No user has asked; additive later. If you choose (b), take the modbases topic's sub-recommendations: no per-code % modified, ignore legacy `Mm`/`Ml`, CHECK on any non-zero inconsistent count (the only flag in the report; conflicts with A1). | MB 1-3, SR D-3 |
| A12 | Initial check thresholds | sign off, or info-only until calibrated | **Moot under A1.** | PR 3 |

## B. Input, performance, build

| # | Question | Options | Recommendation [R] | Sources |
|---|---|---|---|---|
| B1 | Drop protein FASTA? | (a) drop with the tld reader; (b) keep the tld reader for FASTA only | **(a).** Removes `prepare_fasta_fastq_input`, the temp copy and the double gunzip, and enables stdin. "Dropped" means unsupported, not rejected: htslib reads a protein FASTA as nucleotides with non-ACGT letters as N, so the report shows near-100% N. No detection code; the changelog says protein FASTA is unsupported. | PF 1, CL 2, RC 1, plan review |
| B2 | CRAM reference / missing MD | (a) `REF_PATH`/`REF_CACHE` only; (b) add `--reference FILE` (would also let SAMStat compute mismatches for BAMs without MD, e.g. most bwa-mem2/STAR output, and minimap2 without `--MD`, UNVERIFIED) | **(a) for 2.4**, the error message names both variables, and the missing-MD note gives the exact `samtools calmd` / `minimap2 --MD` commands. Keeps "zero new options". **Revisit at G5** with the count of real BAMs lacking MD; if most do, (b) becomes the first 2.5 item. | PF 2, CL 3, SR C15, plan review |
| B3 | libdeflate | (a) keep disabled; (b) auto-detect | **(a).** Builds that differ by machine; gain unmeasured. Revisit with `bench.sh` numbers. | PF 3, SR C25 |
| B4 | Enable bz2/lzma in bundled htslib | (a) yes (in WP5.2, with CRAM); (b) no, "BAM-like CRAMs only" | **(a).** Adds two build deps (also re-added to the bioconda recipe). | PF 4, RC R5 |
| B5 | `-p` and secondary/supplementary counts | (a) count only up to the cutoff; (b) keep 2.3 (depends on batch size) | **(a).** Makes counts independent of batch size and `-t`. In golden regeneration #1. | PF 5 |
| B6 | htslib version, system minimum | (a) bundle latest release (1.24 per research, UNVERIFIED), keep system minimum 1.16; (b) raise minimum | **(a).** Nothing in 2.4 needs a newer API now that modbases is deferred. | RC 2, PF P1, SR C24 |
| B7 | Default `-t` | (a) fixed 4; (b) min(4, cores) | **(a).** Predictable on HPC. `-t` counts all threads: T=4 is 1 reader, 1 decompression thread, 2 workers (formula in WP4.1). | CL, PF P5, SR C13, plan review |
| B8 | CRAM in 2.4 | (a) keep WP5.2, first package dropped on slip; (b) defer to 2.5 "if a user asks" | **(a).** Large cohorts (e.g. nf-core/sarek) deliver CRAM. Cost: 2 build deps, 3 tests, a fixture with reference, reference-error support. | SR §5, plan review |

## C. CLI

| # | Question | Options | Recommendation [R] | Sources |
|---|---|---|---|---|
| C1 | stdin output name | (a) `stdin.samstat.html`; (b) `samstat.html` | **(a).** Keeps the `<name>.samstat.*` glob. | CL 4 |
| C2 | Default output location | (a) beside the input; (b) current directory | **(a)** for 2.x; CWD is a 3.0 question. | CL 5 |
| C3 | Rename `-d` to `-o/--outdir` | (a) keep `-d/--dir`, add `mkdir -p`; (b) rename | **(a).** A rename costs docs and nf-core churn for no capability. | CL, SR C14 |
| C4 | When to remove deprecated options | (a) 3.0; (b) 2.5 | **(a).** Removing an option is breaking. One warning text for all four. | RC R7, CL, SR C16 |

## D. Ecosystem and JSON

| # | Question | Options | Recommendation [R] | Sources |
|---|---|---|---|---|
| D1 | JSON break | (a) clean `format_version` 2, no v1 output, frozen for 2.x; (b) keep v1 and add optional fields | **(a).** Bounded/relative bins cannot be expressed additively. | EC 1, SR C17 |
| D2 | MultiQC error linegraph | (a) in 2.4 (+0.5 d); (b) columns only | **(a).** The only place the differentiator reaches a cohort view. | EC 2, SR §5 |
| D3 | Native MultiQC module | (a) not now; custom content stays; (b) commit to maintaining a Python module | **(a)**, open a module-request issue after bioconda. Your call whether to say in that issue that custom content is permanent. | EC 3 |
| D4 | nf-core module author | (a) you submit and maintain; (b) hand to a volunteer | **Your call.** Recommendation: you submit after 2.4.0 is on bioconda (template `modules/nf-core/cramino`), and ask a volunteer to co-maintain. Not part of the 2.4 code. | EC 4 |
| D5 | Bioconda timing | (a) submit 2.3.0 now, bot bumps to 2.4.0; (b) wait for 2.4.0 | **(a).** If the htslib bump slips, 2.4.0 still ships on 1.16. | EC 5, RC 4, SR C23 |

## E. Testing and release checks

| # | Question | Options | Recommendation [R] | Sources |
|---|---|---|---|---|
| E1 | JSON floats | (a) `%.4f`, from integer sums, no tolerance; (b) 1e-6 tolerance in `compare_golden` | **(a).** Floats appear only in `summary`. Add a tolerance only if a CI-only diff ever appears (`log10` may differ in the last bit between libms, UNVERIFIED). | TS 1, SR C17 |
| E2 | Peak RSS check | (a) manual via `bench.sh`; (b) Linux-only CTest with a ceiling | **(a).** `bench.sh` includes 8 reads of 1 Mb with a 1 GB target; WP4.1 states the memory bound (W+2 batches). A CTest under `/usr/bin/time` or `ulimit -v` is platform-specific and breaks under ASAN. Your call if you want (b). | TS 2, plan review |
| E3 | `dev/` goldens | (a) keep all (428 KB, regenerated twice); (b) cut to a named list | **(a).** Regeneration is one target; they cover samtools edge cases for free. | TS 3 |
| E4 | Real-data check before tagging | (a) required: one real Dorado BAM and one real Revio BAM, by hand; (b) not required | **(a).** Also the precondition for any numeric caption range. **You need to name the two files** (licensing keeps them out of CI). | TS 4, SR D-7 |
| E5 | Size ceilings | JSON/HTML < 1 MB each, or HTML < 2 MB | **JSON < 1 MB, HTML < 1 MB**, HTML tested on the 1 Mb read and on a paired 999 bp file (the short-layout worst case). Holds because there are no per-point tooltips (F11). | TS 5, SR C22, plan review |

## F. Details added by this master plan (please confirm)

| # | Item | Choice made | Why |
|---|---|---|---|
| F1 | N50 and median value above 1 kb | Lower edge of the 1% storage bin; `~` prefix in HTML | Deterministic integer; within 1% as the topics required. Exact below 1000 bp, so a 150 bp library reports 150 (a single 1%-from-100 grid would report 149). |
| F2 | Length plot display bins | One master grid: display edges = exact 1-100, the 1% integer grid up to 996, then exactly the storage edges from 1000. Every display bin is a union of whole storage bins. Used only for plot 1 in long layout and MultiQC; short layout plots storage bins on a linear axis. | Without it, counts jump 10x at the 1 kb storage boundary; a display grid that does not share storage edges gives a sawtooth. |
| F3 | Per-base mismatch split (`mismatches.{A,C,G,T,N}`) | ~~Dropped from JSON v2; one total per bin.~~ **Revised 2026-10-18:** restored as `mismatch_base` (read base), an optional key next to `mismatches`; plot 5 shows one stacked-by-base mismatch row per mapped MAPQ class. | The 2011 paper's headline result (5' added G in CAGE, reaching further into low-MAPQ reads) needs base and MAPQ class together. Read base, not ref->read substitution: it works without MD and is what both papers' examples need; a substitution matrix would roughly double the position block against the 1 MB JSON ceiling (E5). |
| F4 | Two denominators | `aligned` for insertions/deletions, `aligned_md` for mismatches | Keeps mismatch rates correct when some reads lack MD while indel panels always render (SR C10). |
| F5 | Fixed-point scale | `E[q] = round(10^(-q/10) * 2^30)` | No overflow below about 1.7e10 Q0 bases per bin; mean Q above about Q80 is coarse. |
| F6 | Transitional JSON | v2 blocks under `"v2"` in wave 3, promoted in wave 4 | Makes golden regeneration #1 a reviewable, purely additive diff. |
| F7 | Unaligned input | `is_aligned = header has @SQ, or any primary read mapped`; FASTA/FASTQ never. Unaligned: one solid black `all reads` line, no MAPQ table, plots 1-4 only (no N/A lines for 5-6). | Drawing FASTQ or uBAM data as grey "unmapped" would mislead. |
| F8 | Release actions | Tagging, pushing to tld upstream and bioconda PRs are done by you, never by agents | Irreversible and outside the repo. |
| F9 | Workflow size | 23 packages; Opus reviewer only for M/L or Opus-implemented packages (15); 4 gate verifiers: 42 agent runs | Sized by the work, not by a target agent count. |
| F10 | Footer | Version, input basename, citation, secondary/supplementary counts; no command line | A command line leaks absolute paths and differs across `-t`, breaking byte-identity. |
| F11 | Hover and labels | No per-point tooltips, no direct line labels; one series-level `<title>`; MAPQ table is the line key; plot 4 has a 5-swatch legend row | Keeps HTML < 1 MB on short paired input and avoids label-collision code. |
| F12 | Mean Q definition | Summary mean Q and plot 2 over all bases with a quality (clips included, 128-bit internal sum); position plots 3-5 over the aligned part only | One definition per statistic; no uint64 overflow on Tbase runs. |
| F13 | MultiQC position axis | Always the relative (100-bin) axis | Short and long samples in one cohort stay comparable. |
| F14 | Strip slot 4 | `% bases clipped` (aligned only) replaces median length, which stays in the header | Identity ignores clipping; chimeric or adapter-laden reads would otherwise look perfect. One counter, no plot. |

## G. Corrections found while writing this plan

- cli-ux says the `-p` partial note is absent from the HTML. It is present: `stat_report.c:58-59` calls `samstat_partial_warning_section`. 2.4 only moves it into the header line.
- ecosystem §3 writes output names as `<input>.html` / `<input>.json`; the real names are `<input>.samstat.html` and `<input>.samstat.json` (`stat_report.c:263-265`, `docs/multiqc.md`). Unchanged in 2.4.
- presentation and long-reads assume `a->error` lives in a generic struct; it is `struct aln_data` in `src/htsinterface/htsglue.h:24`.
- No current fixture has `=` or X CIGAR operations (checked over `tests/data/*.sam` and `dev/*.sam`), so the WP2.2 X-op change cannot alter existing goldens.
- `v2.3.0` is not tagged yet (`git tag` shows only `v2.0.0`).
