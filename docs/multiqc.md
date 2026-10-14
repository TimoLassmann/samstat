# MultiQC integration

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
or into the `-d` directory:

| file | MultiQC section id | plot type | content |
|---|---|---|---|
| `<file>.samstat_general_mqc.json` | `samstat_general` | `generalstats` | columns in the General Statistics table: reads, % MAPQ >= 30, % unmapped, % properly paired (only if the input has paired reads) |
| `<file>.samstat_mapq_mqc.json` | `samstat_mapq` | `table` | reads per MAPQ bin (`mapq_ge30`, `mapq_10_29`, `mapq_lt10`, `unmapped`), read 1 and read 2 summed |
| `<file>.samstat_length_mqc.json` | `samstat_length` | `linegraph` | read length distribution of read 1 |
| `<file>.samstat_basequal_mqc.json` | `samstat_basequal` | `linegraph` | mean Phred base quality per position (1-based) of read 1 |

Why several files: a MultiQC custom content JSON file holds exactly one section
(one `id`, one `plot_type`). Files from different runs that share an `id` are merged
by MultiQC into one table / plot with one row / line per sample. The sections are
grouped under one `SAMStat` heading (`parent_id: samstat`).

Details:

* Sample name = file name without directory, a trailing `.gz` and one of
  `.sam .bam .fq .fastq .fa .fasta .fna`. It is the key in `data`.
* Percentages are of all counted reads (read 1 + read 2) and written with two decimals.
  Mean qualities have three decimals, computed from the raw quality count arrays over
  all MAPQ bins.
* Read length and quality series longer than 2000 points are summed / averaged in
  blocks of equal width (x is the first value of the block).
* The length and quality files are not written when the input has no read 1 reads.
* Output is deterministic: no timestamps, no absolute paths.
* MAPQ columns use the stable bin ids of the JSON report (`docs/json-format.md`).

## MultiQC sample name cleaning

MultiQC cleans sample names with its own rules (for example a trailing `_qual` is
removed, so `high_qual.sam` appears as `high`). To keep the names SAMStat wrote, add to
`multiqc_config.yaml` (verified with MultiQC 1.35):

```yaml
fn_clean_sample_names: false
```

and run `multiqc -c multiqc_config.yaml <dir>`. Or choose input names MultiQC does not alter.

## Tests

`ctest -L multiqc` checks the files and their ids, plot types and hand-computed values
(`tests/multiqc/`). If `multiqc` is found on `PATH` at configure time, `multiqc:run`
additionally runs MultiQC on three samples and asserts that all appear in the general
statistics and in each SAMStat section. Otherwise that test is not registered
(a STATUS message says so).
