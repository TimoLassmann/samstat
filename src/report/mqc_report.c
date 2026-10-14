#include "tld.h"

#include <inttypes.h>
#include <stdio.h>
#include <string.h>

#include "../plot/plot.h"
#include "../param/param.h"
#include "../collect/collect.h"
#include "../config.h"
#define  MQC_REPORT_IMPORT
#include "mqc_report.h"

/* Layout of the raw counts: see json_report.c. */
#define MQC_N_QUAL 42
#define MQC_N_BIN 4
/* longer series are averaged / summed in blocks so plots stay small */
#define MQC_MAX_POINTS 2000

static const char* const bin_id[MQC_N_BIN] = {"mapq_ge30", "mapq_10_29", "mapq_lt10", "unmapped"};
static const char* const bin_title[MQC_N_BIN] = {"MAPQ >= 30", "MAPQ 10-29", "MAPQ < 10", "Unmapped"};

static int mqc_open(struct samstat_param *p, int id, const char *filename, const char *section, FILE **f);
static int mqc_close(FILE **f);
static int mqc_str(FILE *f, const char *s);
static int mqc_header(FILE *f, const char *id, const char *name, const char *description, const char *plot_type);
static int mqc_sample_name(const char *filename, char **sample);

static int mqc_general(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample);
static int mqc_mapq(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample);
static int mqc_length(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample);
static int mqc_basequal(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample);

int mqc_report_write(struct stat_collection *s, struct samstat_param *p, int id)
{
        char *filename = NULL;
        char *sample = NULL;

        ASSERT(s != NULL, "No stat collection");
        ASSERT(s->mapq_map->n_bin == MQC_N_BIN, "Unexpected number of mapq bins");

        RUN(tld_get_filename(p->infile[id], &filename));
        RUN(mqc_sample_name(filename, &sample));

        RUN(mqc_general(s, p, id, filename, sample));
        RUN(mqc_mapq(s, p, id, filename, sample));
        RUN(mqc_length(s, p, id, filename, sample));
        RUN(mqc_basequal(s, p, id, filename, sample));

        MFREE(sample);
        MFREE(filename);
        return OK;
ERROR:
        MFREE(sample);
        MFREE(filename);
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

        if(len > 3 && strcmp(name + len - 3, ".gz") == 0){
                len -= 3;
                name[len] = 0;
        }
        for(size_t i = 0; i < sizeof(ext) / sizeof(ext[0]); i++){
                size_t l = strlen(ext[i]);
                if(len > l && strcmp(name + len - l, ext[i]) == 0){
                        name[len - l] = 0;
                        break;
                }
        }
        *sample = name;
        return OK;
ERROR:
        return FAIL;
}

int mqc_open(struct samstat_param *p, int id, const char *filename, const char *section, FILE **f)
{
        char path[1024];
        int n;

        if(p->outdir){
                n = snprintf(path, 1024, "%s/%s.samstat_%s_mqc.json", p->outdir, filename, section);
        }else{
                n = snprintf(path, 1024, "%s.samstat_%s_mqc.json", p->infile[id], section);
        }
        ASSERT(n > 0 && n < 1024, "Output path too long");
        LOG_MSG("Writing to: %s", path);
        *f = fopen(path, "w");
        ASSERT(*f != NULL, "Failed to open output file: %s", path);
        return OK;
ERROR:
        return FAIL;
}

int mqc_close(FILE **f)
{
        int r = 0;
        if(*f){
                r = fclose(*f);
                *f = NULL;
        }
        ASSERT(r == 0, "Failed to close MultiQC output file");
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

static void mqc_pct(FILE *f, const char *key, uint64_t num, uint64_t den)
{
        fprintf(f, "\"%s\": %.2f", key, den ? 100.0 * (double) num / (double) den : 0.0);
}

int mqc_general(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample)
{
        FILE *f = NULL;
        uint64_t n_ge30 = 0;
        uint64_t n_unmapped = 0;

        n_ge30 = s->basic_nums[0][0] + (s->n_read2 ? s->basic_nums[1][0] : 0);
        n_unmapped = s->basic_nums[0][3] + (s->n_read2 ? s->basic_nums[1][3] : 0);

        RUN(mqc_open(p, id, filename, "general", &f));
        fprintf(f, "{\n\"id\": \"samstat_general\",\n\"plot_type\": \"generalstats\",\n");
        fprintf(f, "\"pconfig\": [\n");
        fprintf(f, "{\"samstat_total_reads\": {\"title\": \"Reads\", \"description\": \"Total reads counted by SAMStat\", \"namespace\": \"SAMStat\", \"format\": \"{:,.0f}\", \"scale\": \"Blues\"}},\n");
        fprintf(f, "{\"samstat_pct_mapq_ge30\": {\"title\": \"MAPQ >= 30\", \"description\": \"Percent of reads mapped with MAPQ >= 30\", \"namespace\": \"SAMStat\", \"suffix\": \"%%\", \"min\": 0, \"max\": 100, \"format\": \"{:,.1f}\", \"scale\": \"RdYlGn\"}},\n");
        fprintf(f, "{\"samstat_pct_unmapped\": {\"title\": \"Unmapped\", \"description\": \"Percent of reads unmapped\", \"namespace\": \"SAMStat\", \"suffix\": \"%%\", \"min\": 0, \"max\": 100, \"format\": \"{:,.1f}\", \"scale\": \"OrRd\"}}");
        if(s->n_paired){
                fprintf(f, ",\n{\"samstat_pct_proper_pair\": {\"title\": \"Properly paired\", \"description\": \"Percent of reads flagged as properly paired\", \"namespace\": \"SAMStat\", \"suffix\": \"%%\", \"min\": 0, \"max\": 100, \"format\": \"{:,.1f}\", \"scale\": \"RdYlGn\"}}");
        }
        fprintf(f, "\n],\n\"data\": {\n");
        RUN(mqc_str(f, sample));
        fprintf(f, ": {\"samstat_total_reads\": %" PRIu64 ", ", s->n_total);
        mqc_pct(f, "samstat_pct_mapq_ge30", n_ge30, s->n_total);
        fprintf(f, ", ");
        mqc_pct(f, "samstat_pct_unmapped", n_unmapped, s->n_total);
        if(s->n_paired){
                fprintf(f, ", ");
                mqc_pct(f, "samstat_pct_proper_pair", s->n_proper_paired, s->n_total);
        }
        fprintf(f, "}\n}\n}\n");
        RUN(mqc_close(&f));
        return OK;
ERROR:
        if(f){
                fclose(f);
        }
        return FAIL;
}

int mqc_mapq(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample)
{
        FILE *f = NULL;

        RUN(mqc_open(p, id, filename, "mapq", &f));
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
        RUN(mqc_close(&f));
        return OK;
ERROR:
        if(f){
                fclose(f);
        }
        return FAIL;
}

/* Read length distribution of read 1, all mapq bins summed. */
int mqc_length(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample)
{
        FILE *f = NULL;
        struct plot_data *pd = s->len_dist_R1;
        int w = 0;
        int step;
        int first = 1;

        ASSERT(pd->L >= MQC_N_BIN, "Length distribution has too few rows");
        if(s->n_read1 == 0){
                return OK;
        }
        for(int b = 0; b < MQC_N_BIN; b++){
                for(int j = pd->len - 1; j >= w; j--){
                        if(pd->data[b][j]){
                                w = j + 1;
                                break;
                        }
                }
        }
        if(w == 0){
                return OK;
        }
        step = (w + MQC_MAX_POINTS - 1) / MQC_MAX_POINTS;

        RUN(mqc_open(p, id, filename, "length", &f));
        RUN(mqc_header(f, "samstat_length", "Read length distribution",
                       "Number of read 1 reads per read length. Long distributions are summed in blocks; x is the first length of a block.",
                       "linegraph"));
        fprintf(f, "\"pconfig\": {\"id\": \"samstat_length_plot\", \"title\": \"SAMStat: read length distribution\", \"xlab\": \"Read length\", \"ylab\": \"Reads\", \"ymin\": 0, \"xmin\": 0},\n");
        fprintf(f, "\"data\": {\n");
        RUN(mqc_str(f, sample));
        fprintf(f, ": {");
        for(int j = 0; j < w; j += step){
                uint64_t sum = 0;
                for(int k = j; k < j + step && k < w; k++){
                        for(int b = 0; b < MQC_N_BIN; b++){
                                sum += pd->data[b][k];
                        }
                }
                if(step == 1 && sum == 0){
                        continue;
                }
                fprintf(f, "%s\"%d\": %" PRIu64, first ? "" : ", ", j, sum);
                first = 0;
        }
        fprintf(f, "}\n}\n}\n");
        RUN(mqc_close(&f));
        return OK;
ERROR:
        if(f){
                fclose(f);
        }
        return FAIL;
}

/* Mean base quality at each position of read 1 (all mapq bins), from the
   quality count arrays: data[42 * bin + q][position]. */
int mqc_basequal(struct stat_collection *s, struct samstat_param *p, int id, const char *filename, const char *sample)
{
        FILE *f = NULL;
        struct plot_data *pd = s->qual_comp_R1;
        int w = 0;
        int step;
        int first = 1;

        ASSERT(pd->L >= MQC_N_QUAL * MQC_N_BIN, "Quality data has too few rows");
        if(s->n_read1 == 0){
                return OK;
        }
        for(int r = 0; r < MQC_N_QUAL * MQC_N_BIN; r++){
                for(int j = pd->len - 1; j >= w; j--){
                        if(pd->data[r][j]){
                                w = j + 1;
                                break;
                        }
                }
        }
        if(w == 0){
                return OK;
        }
        step = (w + MQC_MAX_POINTS - 1) / MQC_MAX_POINTS;

        RUN(mqc_open(p, id, filename, "basequal", &f));
        RUN(mqc_header(f, "samstat_basequal", "Mean base quality",
                       "Mean Phred base quality per position of read 1. Long reads are averaged in blocks; x is the first position (1-based) of a block.",
                       "linegraph"));
        fprintf(f, "\"pconfig\": {\"id\": \"samstat_basequal_plot\", \"title\": \"SAMStat: mean base quality (read 1)\", \"xlab\": \"Position in read\", \"ylab\": \"Mean Phred quality\", \"ymin\": 0, \"xmin\": 1},\n");
        fprintf(f, "\"data\": {\n");
        RUN(mqc_str(f, sample));
        fprintf(f, ": {");
        for(int j = 0; j < w; j += step){
                uint64_t n = 0;
                uint64_t sum = 0;
                for(int k = j; k < j + step && k < w; k++){
                        for(int b = 0; b < MQC_N_BIN; b++){
                                for(int q = 0; q < MQC_N_QUAL; q++){
                                        uint64_t c = pd->data[MQC_N_QUAL * b + q][k];
                                        n += c;
                                        sum += c * (uint64_t) q;
                                }
                        }
                }
                if(n == 0){
                        continue;
                }
                fprintf(f, "%s\"%d\": %.3f", first ? "" : ", ", j + 1, (double) sum / (double) n);
                first = 0;
        }
        fprintf(f, "}\n}\n}\n");
        RUN(mqc_close(&f));
        return OK;
ERROR:
        if(f){
                fclose(f);
        }
        return FAIL;
}
