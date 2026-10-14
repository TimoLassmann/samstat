# MultiQC integration

> **Set `fn_clean_sample_names: false`.** MultiQC cleans sample names with its own rules
> (for example a trailing `_qual` is removed, so `high_qual.sam` appears as `high`). To keep
> the names SAMStat wrote, add to `multiqc_config.yaml` (verified with MultiQC 1.35):
>
> ```yaml
> fn_clean_sample_names: false
> ```
>
> and run `multiqc -c multiqc_config.yaml <dir>`. Or choose input names MultiQC does not alter.

`samstat --multiqc` writes [MultiQC custom content](https://docs.seqera.io/multiqc/custom_content)
files, so SAMStat results from many runs show up in one MultiQC report.

Validated with MultiQC 1.35 (`uvx multiqc`), three samples, all of them present in the
general statistics table and in each SAMStat section.

## Usage

```
samstat --multiqc -d results/samstat sample1.bam sample2.bam sample3.bam
multiqc results/samstat -o results/multiqc
```

`--multiqc` can be combined with `--json`. Without it no `_mqc` files are written.

## Files

For each input `<file>` (e.g. `sample1.bam`) SAMStat writes, next to the HTML report
or into the `-d` directory (standard input is written as `stdin`):

| file | MultiQC section id | plot type | content |
|---|---|---|---|
| `<file>.samstat_general_mqc.json` | `samstat_general` | `generalstats` | the columns below |
| `<file>.samstat_mapq_mqc.json` | `samstat_mapq` | `table` | reads per MAPQ bin (`mapq_ge30`, `mapq_10_29`, `mapq_lt10`, `unmapped`), read 1 and read 2 summed |
| `<file>.samstat_length_mqc.json` | `samstat_length` | `linegraph` | read 1 reads per read length, logarithmic x axis |
| `<file>.samstat_basequal_mqc.json` | `samstat_basequal` | `linegraph` | read 1 mean Phred quality along the read, all MAPQ classes combined |
| `<file>.samstat_error_mqc.json` | `samstat_error` | `linegraph` | read 1 errors per aligned base along the read, reads with MAPQ >= 30; aligned input only |

### General statistics columns

| id | content |
|---|---|
| `samstat_total_reads` | reads counted (read 1 + read 2) |
| `samstat_pct_mapq_ge30` | percent of reads mapped with MAPQ >= 30 (aligned input only) |
| `samstat_pct_unmapped` | percent of reads unmapped (aligned input only) |
| `samstat_pct_proper_pair` | percent of reads flagged as properly paired (only if the input has paired reads) |
| `samstat_n50` | read length N50 (lower edge of the length bin, exact below 1000 bp) |
| `samstat_mean_q` | mean base quality, over error probabilities (not for FASTA) |
| `samstat_median_identity` | median gap-compressed identity of the aligned reads (only if it can be computed) |

The ids of the first four columns are the ones SAMStat 2.3 wrote and are never renamed.
A column that does not apply to a sample is left out of that sample's file.

Why several files: a MultiQC custom content JSON file holds exactly one section
(one `id`, one `plot_type`). Files from different runs that share an `id` are merged
by MultiQC into one table / plot with one row / line per sample. The sections are
grouped under one `SAMStat` heading (`parent_id: samstat`).

### The two position plots

`samstat_basequal` and `samstat_error` always use the **relative axis**, labelled
`position in read (% of read length)`, for every sample whatever layout its HTML report
has, so short and long samples in one cohort are comparable. The x values are 1 to 100
(a bin ends at that percent of the read; on aligned input the read means its aligned part,
clipped bases excluded). A sample with reads shorter than 100 bases has no data in some
bins and fewer points.

`samstat_error` is `(mismatches per base with MD or X operations) + (insertion and deletion
events per aligned base)`, in percent, for the MAPQ >= 30 reads of read 1. A bin without any
base that carries mismatch information is left out, so a file without `MD` tags and without
`X` operations has no `samstat_error` section. No per-MAPQ plots, composition or identity
histograms are written; those stay in the per-sample HTML report.

### The length plot

x is the lower edge of a display length bin: 1 bp wide up to 100 bp, then bins about 1% wide
(101, 103, ... 986, 996, then the bins of the JSON report from 1000 bp). Every display bin is a
union of whole storage bins of the JSON report, so the counts are exact. Empty bins are left
out, as is the bin of length 0 (the axis is logarithmic).

## Details

* Sample name = file name without directory, a trailing `.gz` and one of
  `.sam .bam .fq .fastq .fa .fasta .fna`. It is the key in `data`.
* Percentages are of all counted reads (read 1 + read 2), written with two decimals;
  median identity has one decimal, mean Q two, mean Q along the read three, error rates four.
  All of them come from the integer counts of the JSON report.
* The length, basequal and error files are not written when there is no read 1 data for them
  (no read 1 reads, no qualities, no aligned bases).
* Output is deterministic: no timestamps, no absolute paths.
* MAPQ columns use the stable bin ids of the JSON report (`docs/json-format.md`).

## Tests

`ctest -L multiqc` checks the files and their ids, plot types and hand-computed values
(`tests/multiqc/`): the exact column list per kind of input, N50 and median identity on
fixtures, the absence of the error section for FASTQ, FASTA and unaligned input, and 100 x
values on the relative axis for a short and a long sample. If `multiqc` is found on `PATH`
at configure time, `multiqc:run` additionally runs MultiQC on three samples and asserts that
all appear in the general statistics and in each SAMStat section, including the error section
(`multiqc:run_error_section`). Otherwise those tests are not registered (a STATUS message says so).
