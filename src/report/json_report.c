#include "tld.h"

#include <inttypes.h>
#include <stdio.h>
#include <string.h>

#include "../plot/plot.h"
#include "../param/param.h"
#include "../collect/collect.h"
#include "../config.h"
#define  JSON_REPORT_IMPORT
#include "json_report.h"

/* Layout of the raw counts in struct stat_collection:
 *   len_dist      data[bin][length]
 *   base_comp     data[5  * bin + base][position]   (A,C,G,T,N)
 *   qual_comp     data[42 * bin + q][position]
 *   mis           data[5  * bin + base][position]
 *   ins / del     data[bin][position]
 */
#define JSON_N_BASE 5
#define JSON_N_QUAL 42
#define JSON_N_BIN 4

static const char* const bin_id[JSON_N_BIN] = {"mapq_ge30", "mapq_10_29", "mapq_lt10", "unmapped"};
static const char* const base_name[JSON_N_BASE] = {"A", "C", "G", "T", "N"};

/* Output is streamed: the buffer is flushed to the file whenever it holds
   JSON_FLUSH_BYTES, so memory stays small for very long reads. */
#define JSON_FLUSH_BYTES 32768
#define JSON_BUF_BYTES 65536
struct jout {
        tld_strbuf *b;
        FILE *f;
};

static int json_append(struct jout *o, char *s);
static int json_char(struct jout *o, char c);
static int json_flush(struct jout *o);
static int json_str(struct jout *o, const char *s);
static int json_u64(struct jout *o, uint64_t v);
static int json_key(struct jout *o, const char *key);
static int json_array(struct jout *o, uint64_t *row, int n);
static int json_width(uint64_t **rows, int n_rows, int len);

static int json_bin_counts(struct jout *o, uint64_t *counts);
static int json_read(struct jout *o, struct stat_collection *s, int read);
static int json_section_simple(struct jout *o, struct plot_data *pd, const char *name);
static int json_section_grouped(struct jout *o, struct plot_data *pd, const char *name, int group, int is_qual);

int json_report_write(struct stat_collection *s, struct samstat_param *p, int id)
{
        struct jout out = {NULL, NULL};
        struct jout *o = &out;
        char *filename = NULL;
        char path[1024];
        int n;

        ASSERT(s != NULL, "No stat collection");
        ASSERT(s->mapq_map->n_bin == JSON_N_BIN, "Unexpected number of mapq bins");

        RUN(tld_get_filename(p->infile[id], &filename));

        if(p->outdir){
                n = snprintf(path, 1024, "%s/%s.samstat.json", p->outdir, filename);
        }else{
                n = snprintf(path, 1024, "%s.samstat.json", p->infile[id]);
        }
        ASSERT(n > 0 && n < 1024, "Output path too long");
        LOG_MSG("Writing to: %s", path);
        out.f = fopen(path, "w");
        ASSERT(out.f != NULL, "Failed to open output file: %s", path);
        RUN(tld_strbuf_alloc(&out.b, JSON_BUF_BYTES));

        RUN(json_append(o, "{\n\"format_version\": 1,\n"));
        RUN(json_append(o, "\"meta\": {\"samstat_version\": "));
        RUN(json_str(o, PACKAGE_VERSION));
        RUN(json_append(o, ", \"input_file\": "));
        RUN(json_str(o, filename));
        RUN(json_append(o, "},\n"));

        RUN(json_append(o, "\"summary\": {\n"));
        RUN(json_append(o, "\"total_reads\": "));
        RUN(json_u64(o, s->n_total));
        RUN(json_append(o, ",\n\"read1_reads\": "));
        RUN(json_u64(o, s->n_read1));
        RUN(json_append(o, ",\n\"read2_reads\": "));
        RUN(json_u64(o, s->n_read2));
        RUN(json_append(o, ",\n\"paired\": "));
        RUN(json_u64(o, s->n_paired));
        RUN(json_append(o, ",\n\"properly_paired\": "));
        RUN(json_u64(o, s->n_proper_paired));
        RUN(json_append(o, ",\n\"secondary\": "));
        RUN(json_u64(o, s->n_secondary));
        RUN(json_append(o, ",\n\"supplementary\": "));
        RUN(json_u64(o, s->n_supplementary));
        RUN(json_append(o, ",\n\"is_aligned\": "));
        RUN(json_append(o, s->is_aligned ? "true" : "false"));
        RUN(json_append(o, ",\n\"is_partial_report\": "));
        RUN(json_append(o, s->is_partial_report ? "true" : "false"));
        RUN(json_append(o, ",\n\"mapq_counts\": {\n\"read1\": "));
        RUN(json_bin_counts(o, s->basic_nums[0]));
        if(s->n_read2){
                RUN(json_append(o, ",\n\"read2\": "));
                RUN(json_bin_counts(o, s->basic_nums[1]));
        }
        RUN(json_append(o, "\n}\n},\n"));

        RUN(json_append(o, "\"mapq_bins\": [\n"));
        for(int i = 0; i < JSON_N_BIN; i++){
                RUN(json_append(o, "{\"id\": "));
                RUN(json_str(o, bin_id[i]));
                RUN(json_append(o, ", \"label\": "));
                RUN(json_str(o, s->mapq_map->description[i]));
                RUN(json_append(o, i + 1 < JSON_N_BIN ? "},\n" : "}\n"));
        }
        RUN(json_append(o, "],\n"));

        RUN(json_append(o, "\"read1\": "));
        RUN(json_read(o, s, 1));
        if(s->n_read2){
                RUN(json_append(o, ",\n\"read2\": "));
                RUN(json_read(o, s, 2));
        }
        RUN(json_append(o, "\n}\n"));

        RUN(json_flush(o));
        ASSERT(fclose(out.f) == 0, "Failed to close %s", path);
        out.f = NULL;

        tld_strbuf_free(out.b);
        MFREE(filename);
        return OK;
ERROR:
        if(out.f){
                fclose(out.f);
        }
        tld_strbuf_free(out.b);
        MFREE(filename);
        return FAIL;
}

int json_append(struct jout *o, char *s)
{
        RUN(tld_append(o->b, s));
        if(o->b->len >= JSON_FLUSH_BYTES){
                RUN(json_flush(o));
        }
        return OK;
ERROR:
        return FAIL;
}

int json_char(struct jout *o, char c)
{
        RUN(tld_append_char(o->b, c));
        if(o->b->len >= JSON_FLUSH_BYTES){
                RUN(json_flush(o));
        }
        return OK;
ERROR:
        return FAIL;
}

int json_flush(struct jout *o)
{
        size_t len = (size_t) o->b->len;
        ASSERT(fwrite(TLD_STR(o->b), 1, len, o->f) == len, "Failed to write JSON report");
        RUN(tld_strbuf_clear(o->b));
        return OK;
ERROR:
        return FAIL;
}

int json_read(struct jout *o, struct stat_collection *s, int read)
{
        struct plot_data *len_d = s->len_dist_R1;
        struct plot_data *base_d = s->base_comp_R1;
        struct plot_data *qual_d = s->qual_comp_R1;
        struct plot_data *mis_d = s->mis_R1;
        struct plot_data *ins_d = s->ins_R1;
        struct plot_data *del_d = s->del_R1;
        uint64_t n = s->n_read1;

        if(read == 2){
                len_d = s->len_dist_R2;
                base_d = s->base_comp_R2;
                qual_d = s->qual_comp_R2;
                mis_d = s->mis_R2;
                ins_d = s->ins_R2;
                del_d = s->del_R2;
                n = s->n_read2;
        }

        RUN(json_append(o, "{\n\"n_reads\": "));
        RUN(json_u64(o, n));
        RUN(json_append(o, ",\n"));
        RUN(json_section_simple(o, len_d, "length"));
        RUN(json_append(o, ",\n"));
        RUN(json_section_grouped(o, base_d, "base_composition", JSON_N_BASE, 0));
        RUN(json_append(o, ",\n"));
        RUN(json_section_grouped(o, qual_d, "base_quality", JSON_N_QUAL, 1));
        RUN(json_append(o, ",\n"));
        RUN(json_section_grouped(o, mis_d, "mismatches", JSON_N_BASE, 0));
        RUN(json_append(o, ",\n"));
        RUN(json_section_simple(o, ins_d, "insertions"));
        RUN(json_append(o, ",\n"));
        RUN(json_section_simple(o, del_d, "deletions"));
        RUN(json_append(o, "\n}"));
        return OK;
ERROR:
        return FAIL;
}

/* One array per mapq bin: pd->data[bin][i], trimmed after the last non-zero entry. */
int json_section_simple(struct jout *o, struct plot_data *pd, const char *name)
{
        ASSERT(pd->L >= JSON_N_BIN, "Plot data %s has too few rows", name);
        RUN(json_key(o, name));
        RUN(json_append(o, "{"));
        for(int b = 0; b < JSON_N_BIN; b++){
                int w = json_width(pd->data + b, 1, pd->len);
                RUN(json_key(o, bin_id[b]));
                RUN(json_array(o, pd->data[b], w));
                if(b + 1 < JSON_N_BIN){
                        RUN(json_append(o, ", "));
                }
        }
        RUN(json_append(o, "}"));
        return OK;
ERROR:
        return FAIL;
}

/* Per bin: 'group' rows (bases or Q values) of pd->data[group * bin + k][position].
   Bases: object keyed by base. Quality: array indexed by Q value. All rows of a
   bin share one width (last non-zero position over all rows, plus one). */
int json_section_grouped(struct jout *o, struct plot_data *pd, const char *name, int group, int is_qual)
{
        ASSERT(pd->L >= group * JSON_N_BIN, "Plot data %s has too few rows", name);
        RUN(json_key(o, name));
        RUN(json_append(o, "{\n"));
        for(int b = 0; b < JSON_N_BIN; b++){
                uint64_t **rows = pd->data + group * b;
                int w = json_width(rows, group, pd->len);
                RUN(json_key(o, bin_id[b]));
                if(is_qual){
                        if(w == 0){
                                RUN(json_append(o, "[]"));
                        }else{
                                RUN(json_append(o, "["));
                                for(int k = 0; k < group; k++){
                                        RUN(json_array(o, rows[k], w));
                                        if(k + 1 < group){
                                                RUN(json_append(o, ", "));
                                        }
                                }
                                RUN(json_append(o, "]"));
                        }
                }else{
                        RUN(json_append(o, "{"));
                        for(int k = 0; k < group; k++){
                                RUN(json_key(o, base_name[k]));
                                RUN(json_array(o, rows[k], w));
                                if(k + 1 < group){
                                        RUN(json_append(o, ", "));
                                }
                        }
                        RUN(json_append(o, "}"));
                }
                RUN(json_append(o, b + 1 < JSON_N_BIN ? ",\n" : "\n"));
        }
        RUN(json_append(o, "}"));
        return OK;
ERROR:
        return FAIL;
}

int json_bin_counts(struct jout *o, uint64_t *counts)
{
        RUN(json_append(o, "{"));
        for(int b = 0; b < JSON_N_BIN; b++){
                RUN(json_key(o, bin_id[b]));
                RUN(json_u64(o, counts[b]));
                if(b + 1 < JSON_N_BIN){
                        RUN(json_append(o, ", "));
                }
        }
        RUN(json_append(o, "}"));
        return OK;
ERROR:
        return FAIL;
}

/* number of leading entries needed to cover every non-zero count */
int json_width(uint64_t **rows, int n_rows, int len)
{
        int w = 0;
        for(int i = 0; i < n_rows; i++){
                for(int j = len - 1; j >= w; j--){
                        if(rows[i][j]){
                                w = j + 1;
                                break;
                        }
                }
        }
        return w;
}

int json_array(struct jout *o, uint64_t *row, int n)
{
        RUN(json_append(o, "["));
        for(int i = 0; i < n; i++){
                if(i){
                        RUN(json_append(o, ","));
                }
                RUN(json_u64(o, row[i]));
        }
        RUN(json_append(o, "]"));
        return OK;
ERROR:
        return FAIL;
}

int json_u64(struct jout *o, uint64_t v)
{
        char buf[32];
        snprintf(buf, 32, "%" PRIu64, v);
        RUN(json_append(o, buf));
        return OK;
ERROR:
        return FAIL;
}

int json_key(struct jout *o, const char *key)
{
        RUN(json_str(o, key));
        RUN(json_append(o, ": "));
        return OK;
ERROR:
        return FAIL;
}

int json_str(struct jout *o, const char *s)
{
        char buf[8];
        RUN(json_char(o, '"'));
        for(; *s; s++){
                unsigned char c = (unsigned char) *s;
                switch(c){
                case '"':
                        RUN(json_append(o, "\\\""));
                        break;
                case '\\':
                        RUN(json_append(o, "\\\\"));
                        break;
                case '\n':
                        RUN(json_append(o, "\\n"));
                        break;
                case '\r':
                        RUN(json_append(o, "\\r"));
                        break;
                case '\t':
                        RUN(json_append(o, "\\t"));
                        break;
                default:
                        if(c < 0x20){
                                snprintf(buf, 8, "\\u%04x", c);
                                RUN(json_append(o, buf));
                        }else{
                                RUN(json_char(o, (char) c));
                        }
                }
        }
        RUN(json_char(o, '"'));
        return OK;
ERROR:
        return FAIL;
}
