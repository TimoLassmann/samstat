#include "tldlite.h"

#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include "../param/param.h"
#include "../collect/collect.h"
#include "config.h"
#define  MQC_REPORT_IMPORT
#include "mqc_report.h"
#include "outfile.h"

#define MQC_N_BIN 4

static const char* const bin_id[MQC_N_BIN] = {"mapq_ge30", "mapq_10_29", "mapq_lt10", "unmapped"};
static const char* const bin_title[MQC_N_BIN] = {"MAPQ >= 30", "MAPQ 10-29", "MAPQ < 10", "Unmapped"};

static int mqc_open(struct samstat_param *p, int id, const char *filename, const char *section, struct outfile **o);
static int mqc_close(struct outfile **o);
static int mqc_str(FILE *f, const char *s);
static int mqc_header(FILE *f, const char *id, const char *name, const char *description, const char *plot_type);
static int mqc_sample_name(const char *filename, char **sample);

static int mqc_general(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample);
static int mqc_mapq(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample);
static int mqc_length(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample);
static int mqc_basequal(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample);
static int mqc_error(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample);

int mqc_report_write(struct stat_collection *s, struct samstat_param *p, int id)
{
        char *filename = NULL;
        char *sample = NULL;

        ASSERT(s != NULL, "No stat collection");
        ASSERT(s->mapq_map->n_bin == MQC_N_BIN, "Unexpected number of mapq bins");

        RUN(tld_get_filename((char *) samstat_output_stem(p, id), &filename));
        RUN(mqc_sample_name(filename, &sample));

        RUN(mqc_general(s, p, id, filename, sample));
        RUN(mqc_mapq(s, p, id, filename, sample));
        RUN(mqc_length(s, p, id, filename, sample));
        RUN(mqc_basequal(s, p, id, filename, sample));
        RUN(mqc_error(s, p, id, filename, sample));

        MFREE(sample);
        MFREE(filename);
        return OK;
ERROR:
        if(sample){
                MFREE(sample);
        }
        if(filename){
                MFREE(filename);
        }
        return FAIL;
}

/* Sample name = file name without directory and without a trailing
   .gz and one known sequence / alignment extension. */
int mqc_sample_name(const char *filename, char **sample)
{
        static const char* const ext[] = {".sam", ".bam", ".cram", ".fq", ".fastq", ".fa", ".fasta", ".fna"};
        char *name = NULL;
        size_t len = strlen(filename);

        MMALLOC(name, len + 1);
        memcpy(name, filename, len + 1);

        if(len > 3 && strcasecmp(name + len - 3, ".gz") == 0){
                len -= 3;
                name[len] = 0;
        }
        for(size_t i = 0; i < sizeof(ext) / sizeof(ext[0]); i++){
                size_t l = strlen(ext[i]);
                if(len > l && strcasecmp(name + len - l, ext[i]) == 0){
                        name[len - l] = 0;
                        break;
                }
        }
        *sample = name;
        return OK;
ERROR:
        return FAIL;
}

/* Each file is written to a temporary file and renamed when complete (outfile.h):
   a failed write leaves no truncated _mqc file for MultiQC to pick up. */
int mqc_open(struct samstat_param *p, int id, const char *filename, const char *section, struct outfile **o)
{
        char *path = NULL;

        if(p->outdir){
                RUN(outfile_path(&path, "%s/%s.samstat_%s_mqc.json", p->outdir, filename, section));
        }else{
                RUN(outfile_path(&path, "%s.samstat_%s_mqc.json", samstat_output_stem(p, id), section));
        }
        if(p->verbose){
                LOG_MSG("Writing to: %s", path);
        }
        RUN(outfile_open(o, path));
        MFREE(path);
        return OK;
ERROR:
        if(path){
                MFREE(path);
        }
        return FAIL;
}

/* The fprintf/fputc calls of the sections are not checked one by one: the error
   flag of the stream is sticky and outfile_commit checks it before the rename. */
int mqc_close(struct outfile **o)
{
        RUN(outfile_commit(o));
        return OK;
ERROR:
        return FAIL;
}

int mqc_str(FILE *f, const char *s)
{
        fputc('"', f);
        for(; *s; s++){
                unsigned char c = (unsigned char) *s;
                switch(c){
                case '"':
                        fputs("\\\"", f);
                        break;
                case '\\':
                        fputs("\\\\", f);
                        break;
                default:
                        if(c < 0x20){
                                fprintf(f, "\\u%04x", c);
                        }else{
                                fputc(c, f);
                        }
                }
        }
        fputc('"', f);
        return OK;
}

/* Common leading keys. Sections share a parent so MultiQC groups them under 'SAMStat'. */
int mqc_header(FILE *f, const char *id, const char *name, const char *description, const char *plot_type)
{
        fprintf(f, "{\n\"id\": ");
        RUN(mqc_str(f, id));
        fprintf(f, ",\n\"parent_id\": \"samstat\",\n\"parent_name\": \"SAMStat\",\n\"parent_description\": \"Statistics of SAM/BAM/FASTQ files.\"");
        fprintf(f, ",\n\"section_name\": ");
        RUN(mqc_str(f, name));
        fprintf(f, ",\n\"description\": ");
        RUN(mqc_str(f, description));
        fprintf(f, ",\n\"plot_type\": ");
        RUN(mqc_str(f, plot_type));
        fprintf(f, ",\n");
        return OK;
ERROR:
        return FAIL;
}

/* General statistics (standard-plots.md section 8). A column that does not apply to the
   input is left out of both pconfig and data: pct_mapq_ge30 and pct_unmapped for input
   that is not aligned, pct_proper_pair for unpaired input, n50 without sequence, mean_q without
   qualities, median_identity without a computable identity. */
int mqc_general(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample)
{
        struct outfile *of = NULL;
        FILE *f = NULL;
        int has_ge30 = !COLLECT_IS_NA(s->pct_mapq_ge30);
        int has_unmapped = !COLLECT_IS_NA(s->pct_mapped);
        int has_pair = s->n_paired > 0;
        int has_n50 = s->n_total > 0 && s->total_bases > 0;
        int has_q = !COLLECT_IS_NA(s->mean_q);
        int has_ident = !COLLECT_IS_NA(s->median_identity);

        RUN(mqc_open(p, id, filename, "general", &of));
        f = of->f;
        fprintf(f, "{\n\"id\": \"samstat_general\",\n\"plot_type\": \"generalstats\",\n");
        fprintf(f, "\"pconfig\": [\n");
        fprintf(f, "{\"samstat_total_reads\": {\"title\": \"Reads\", \"description\": \"Total reads counted by SAMStat\", \"namespace\": \"SAMStat\", \"format\": \"{:,.0f}\", \"scale\": \"Blues\"}}");
        if(has_ge30){
                fprintf(f, ",\n{\"samstat_pct_mapq_ge30\": {\"title\": \"MAPQ >= 30\", \"description\": \"Percent of reads mapped with MAPQ >= 30\", \"namespace\": \"SAMStat\", \"suffix\": \"%%\", \"min\": 0, \"max\": 100, \"format\": \"{:,.1f}\", \"scale\": \"RdYlGn\"}}");
        }
        if(has_unmapped){
                fprintf(f, ",\n{\"samstat_pct_unmapped\": {\"title\": \"Unmapped\", \"description\": \"Percent of reads unmapped\", \"namespace\": \"SAMStat\", \"suffix\": \"%%\", \"min\": 0, \"max\": 100, \"format\": \"{:,.1f}\", \"scale\": \"OrRd\"}}");
        }
        if(has_pair){
                fprintf(f, ",\n{\"samstat_pct_proper_pair\": {\"title\": \"Properly paired\", \"description\": \"Percent of reads flagged as properly paired\", \"namespace\": \"SAMStat\", \"suffix\": \"%%\", \"min\": 0, \"max\": 100, \"format\": \"{:,.1f}\", \"scale\": \"RdYlGn\"}}");
        }
        if(has_n50){
                fprintf(f, ",\n{\"samstat_n50\": {\"title\": \"N50\", \"description\": \"Read length N50 (bp; exact below 1 kb, otherwise within 1%%)\", \"namespace\": \"SAMStat\", \"suffix\": \" bp\", \"min\": 0, \"format\": \"{:,.0f}\", \"scale\": \"Greens\"}}");
        }
        if(has_q){
                fprintf(f, ",\n{\"samstat_mean_q\": {\"title\": \"Mean Q\", \"description\": \"Mean base quality (Phred, from the mean error probability)\", \"namespace\": \"SAMStat\", \"min\": 0, \"format\": \"{:,.1f}\", \"scale\": \"RdYlGn\"}}");
        }
        if(has_ident){
                fprintf(f, ",\n{\"samstat_median_identity\": {\"title\": \"Median identity\", \"description\": \"Median gap-compressed identity of aligned reads\", \"namespace\": \"SAMStat\", \"suffix\": \"%%\", \"min\": 0, \"max\": 100, \"format\": \"{:,.1f}\", \"scale\": \"RdYlGn\"}}");
        }
        fprintf(f, "\n],\n\"data\": {\n");
        RUN(mqc_str(f, sample));
        fprintf(f, ": {\"samstat_total_reads\": %" PRIu64, s->n_total);
        if(has_ge30){
                fprintf(f, ", \"samstat_pct_mapq_ge30\": %.2f", s->pct_mapq_ge30);
        }
        if(has_unmapped){
                fprintf(f, ", \"samstat_pct_unmapped\": %.2f", 100.0 - s->pct_mapped);
        }
        if(has_pair){
                fprintf(f, ", \"samstat_pct_proper_pair\": %.2f", s->n_total ? 100.0 * (double) s->n_proper_paired / (double) s->n_total : 0.0);
        }
        if(has_n50){
                fprintf(f, ", \"samstat_n50\": %" PRIu64, s->n50_best);
        }
        if(has_q){
                fprintf(f, ", \"samstat_mean_q\": %.2f", s->mean_q);
        }
        if(has_ident){
                fprintf(f, ", \"samstat_median_identity\": %.1f", s->median_identity);
        }
        fprintf(f, "}\n}\n}\n");
        RUN(mqc_close(&of));
        return OK;
ERROR:
        outfile_abort(&of);
        return FAIL;
}

int mqc_mapq(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample)
{
        struct outfile *of = NULL;
        FILE *f = NULL;

        RUN(mqc_open(p, id, filename, "mapq", &of));
        f = of->f;
        RUN(mqc_header(f, "samstat_mapq", "Mapping quality",
                       "Number of reads per mapping quality bin (read 1 and read 2 combined). Unmapped also holds all FASTA/FASTQ reads.",
                       "table"));
        fprintf(f, "\"pconfig\": {\"id\": \"samstat_mapq_table\", \"title\": \"SAMStat: reads per MAPQ bin\", \"sort_rows\": false},\n");
        fprintf(f, "\"headers\": {\n");
        for(int b = 0; b < MQC_N_BIN; b++){
                fprintf(f, "\"%s\": {\"title\": \"%s\", \"description\": \"Reads in bin: %s\", \"format\": \"{:,.0f}\", \"scale\": \"Blues\"}%s\n",
                        bin_id[b], bin_title[b], bin_title[b], b + 1 < MQC_N_BIN ? "," : "");
        }
        fprintf(f, "},\n\"data\": {\n");
        RUN(mqc_str(f, sample));
        fprintf(f, ": {");
        for(int b = 0; b < MQC_N_BIN; b++){
                uint64_t v = s->basic_nums[0][b] + (s->n_read2 ? s->basic_nums[1][b] : 0);
                fprintf(f, "\"%s\": %" PRIu64 "%s", bin_id[b], v, b + 1 < MQC_N_BIN ? ", " : "");
        }
        fprintf(f, "}\n}\n}\n");
        RUN(mqc_close(&of));
        return OK;
ERROR:
        outfile_abort(&of);
        return FAIL;
}

/* Read length distribution of read 1, all MAPQ classes summed, on the display bins of
   standard-plots.md 3.1 (every display bin is a union of whole storage bins, so the sums
   are exact). x = lower edge of the display bin; empty bins are left out; the bin of
   length 0 is left out because the axis is logarithmic. */
int mqc_length(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample)
{
        struct outfile *of = NULL;
        FILE *f = NULL;
        uint64_t *sum = NULL;
        int n_disp = collect_n_display_bins();
        int first = 1;

        if(s->n_read1 == 0){
                return OK;
        }
        MMALLOC(sum, sizeof(uint64_t) * (size_t) n_disp);
        for(int i = 0; i < n_disp; i++){
                sum[i] = 0;
        }
        for(int b = 0; b < COLLECT_N_LEN_BINS; b++){
                int d = collect_display_bin(b);
                for(int m = 0; m < COLLECT_N_MAPQ; m++){
                        sum[d] += s->acc->len_reads[0][m][b];
                }
        }

        RUN(mqc_open(p, id, filename, "length", &of));
        f = of->f;
        RUN(mqc_header(f, "samstat_length", "Read length distribution",
                       "Number of read 1 reads per read length. Lengths above 100 bp are grouped into bins about 1% wide; x is the lower edge of a bin.",
                       "linegraph"));
        fprintf(f, "\"pconfig\": {\"id\": \"samstat_length_plot\", \"title\": \"SAMStat: read length distribution\", \"xlab\": \"Read length (bp)\", \"ylab\": \"Reads\", \"ymin\": 0, \"xLog\": true},\n");
        fprintf(f, "\"data\": {\n");
        RUN(mqc_str(f, sample));
        fprintf(f, ": {");
        for(int i = 1; i < n_disp; i++){
                if(sum[i]){
                        fprintf(f, "%s\"%" PRIu64 "\": %" PRIu64, first ? "" : ", ", collect_display_edge(i), sum[i]);
                        first = 0;
                }
        }
        fprintf(f, "}\n}\n}\n");
        RUN(mqc_close(&of));
        MFREE(sum);
        return OK;
ERROR:
        outfile_abort(&of);
        if(sum){
                MFREE(sum);
        }
        return FAIL;
}

#define MQC_X_LABEL "position in read (%% of read length)"

/* Mean Q of read 1 on the relative position axis (100 bins), all MAPQ classes combined,
   from the integer error sums. x = 1..100 (the bin ends at that percent of the read). */
int mqc_basequal(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample)
{
        struct outfile *of = NULL;
        FILE *f = NULL;
        int first = 1;
        int any = 0;

        if(s->n_read1 == 0){
                return OK;
        }
        for(int i = 0; i < COLLECT_REL_BINS && !any; i++){
                for(int m = 0; m < COLLECT_N_MAPQ; m++){
                        if(s->acc->pos_rel[0][m][i].qual_bases){
                                any = 1;
                        }
                }
        }
        if(!any){
                return OK;
        }

        RUN(mqc_open(p, id, filename, "basequal", &of));
        f = of->f;
        RUN(mqc_header(f, "samstat_basequal", "Mean base quality",
                       "Mean Phred base quality of read 1 along the read, in 100 bins of the read length (aligned part for aligned input). The same axis for every sample.",
                       "linegraph"));
        fprintf(f, "\"pconfig\": {\"id\": \"samstat_basequal_plot\", \"title\": \"SAMStat: mean base quality (read 1)\", \"xlab\": \"" MQC_X_LABEL "\", \"ylab\": \"Mean Phred quality\", \"ymin\": 0, \"xmin\": 1, \"xmax\": 100},\n");
        fprintf(f, "\"data\": {\n");
        RUN(mqc_str(f, sample));
        fprintf(f, ": {");
        for(int i = 0; i < COLLECT_REL_BINS; i++){
                collect_u128 err = 0;
                uint64_t n = 0;
                double q;
                for(int m = 0; m < COLLECT_N_MAPQ; m++){
                        err += s->acc->pos_rel[0][m][i].qual_error_sum;
                        n += s->acc->pos_rel[0][m][i].qual_bases;
                }
                q = collect_mean_q(err, n);
                if(COLLECT_IS_NA(q)){
                        continue;
                }
                fprintf(f, "%s\"%d\": %.3f", first ? "" : ", ", i + 1, q);
                first = 0;
        }
        fprintf(f, "}\n}\n}\n");
        RUN(mqc_close(&of));
        return OK;
ERROR:
        outfile_abort(&of);
        return FAIL;
}

/* Error rate of read 1 reads with MAPQ >= 30 on the relative axis: (mismatches per base
   with MD or X) + (insertion and deletion events per aligned base), in percent. A bin
   without any base that carries mismatch information is left out, as is the whole file
   for input that is not aligned or has no such bin. */
int mqc_error(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample)
{
        struct outfile *of = NULL;
        FILE *f = NULL;
        const struct collect_pos_bin *bins;
        int first = 1;
        int any = 0;

        if(!s->is_aligned || s->n_read1 == 0){
                return OK;
        }
        bins = s->acc->pos_rel[0][COLLECT_MAPQ_GE30];
        for(int i = 0; i < COLLECT_REL_BINS; i++){
                if(bins[i].aligned && bins[i].aligned_md){
                        any = 1;
                        break;
                }
        }
        if(!any){
                return OK;
        }

        RUN(mqc_open(p, id, filename, "error", &of));
        f = of->f;
        RUN(mqc_header(f, "samstat_error", "Error rate along the read",
                       "Mismatches, insertions and deletions per aligned base of read 1 reads with MAPQ >= 30, in 100 bins of the read length (aligned part). Mismatches come from MD or CIGAR X; each gap opening counts once.",
                       "linegraph"));
        fprintf(f, "\"pconfig\": {\"id\": \"samstat_error_plot\", \"title\": \"SAMStat: errors along the read (read 1, MAPQ >= 30)\", \"xlab\": \"" MQC_X_LABEL "\", \"ylab\": \"Errors per aligned base (%%)\", \"ymin\": 0, \"xmin\": 1, \"xmax\": 100},\n");
        fprintf(f, "\"data\": {\n");
        RUN(mqc_str(f, sample));
        fprintf(f, ": {");
        for(int i = 0; i < COLLECT_REL_BINS; i++){
                double e;
                if(!bins[i].aligned || !bins[i].aligned_md){
                        continue;
                }
                e = 100.0 * ((double) bins[i].mismatches / (double) bins[i].aligned_md
                             + (double) (bins[i].insertions + bins[i].deletions) / (double) bins[i].aligned);
                fprintf(f, "%s\"%d\": %.4f", first ? "" : ", ", i + 1, e);
                first = 0;
        }
        fprintf(f, "}\n}\n}\n");
        RUN(mqc_close(&of));
        return OK;
ERROR:
        outfile_abort(&of);
        return FAIL;
}
