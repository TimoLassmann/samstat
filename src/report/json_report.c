#include "tldlite.h"

#include <inttypes.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../param/param.h"
#include "../collect/collect.h"
#include "config.h"
#define  JSON_REPORT_IMPORT
#include "json_report.h"
#include "outfile.h"

#define JSON_N_BASE 5
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

static int json_bin_counts(struct jout *o, uint64_t *counts);
static int json_double(struct jout *o, double v);
static int json_strided(struct jout *o, const char *base, size_t stride, int n);
static int json_v2(struct jout *o, struct stat_collection *s, const char *filename);
static int json_v2_summary(struct jout *o, struct stat_collection *s);
static int json_v2_read(struct jout *o, struct stat_collection *s, int r);
static int json_v2_per_class(struct jout *o, const char *base, size_t class_stride, size_t stride, int n);
static int json_v2_per_class_base(struct jout *o, const char *base, size_t class_stride, size_t stride, int n);
static int json_v2_pos_block(struct jout *o, const struct collect_pos_bin *bins, int n_bins);

int json_report_write(struct stat_collection *s, struct samstat_param *p, int id)
{
        struct jout out = {NULL, NULL};
        struct jout *o = &out;
        struct outfile *of = NULL;
        char *filename = NULL;
        char *path = NULL;

        ASSERT(s != NULL, "No stat collection");
        ASSERT(s->mapq_map->n_bin == JSON_N_BIN, "Unexpected number of mapq bins");

        RUN(tld_get_filename((char *) samstat_output_stem(p, id), &filename));

        if(p->outdir){
                RUN(outfile_path(&path, "%s/%s.samstat.json", p->outdir, filename));
        }else{
                RUN(outfile_path(&path, "%s.samstat.json", samstat_output_stem(p, id)));
        }
        if(p->verbose){
                LOG_MSG("Writing to: %s", path);
        }
        RUN(outfile_open(&of, path));
        out.f = of->f;
        RUN(tld_strbuf_alloc(&out.b, JSON_BUF_BYTES));

        RUN(json_v2(o, s, filename));
        RUN(json_append(o, "\n"));

        RUN(json_flush(o));
        RUN(outfile_commit(&of));

        tld_strbuf_free(out.b);
        MFREE(path);
        MFREE(filename);
        return OK;
ERROR:
        outfile_abort(&of);
        tld_strbuf_free(out.b);
        if(path){
                MFREE(path);
        }
        if(filename){
                MFREE(filename);
        }
        return FAIL;
}

/* ---- format_version 2 (standard-plots.md section 7) ---- */

int json_v2(struct jout *o, struct stat_collection *s, const char *filename)
{
        ASSERT(s->acc != NULL, "No v2 accumulators");
        RUN(json_append(o, "{\n\"format_version\": 2,\n"));
        RUN(json_append(o, "\"meta\": {\"samstat_version\": "));
        RUN(json_str(o, PACKAGE_VERSION));
        RUN(json_append(o, ", \"input_file\": "));
        RUN(json_str(o, filename));
        RUN(json_append(o, "},\n"));
        RUN(json_v2_summary(o, s));
        RUN(json_append(o, ",\n\"mapq_bins\": [\n"));
        for(int i = 0; i < JSON_N_BIN; i++){
                RUN(json_append(o, "{\"id\": "));
                RUN(json_str(o, bin_id[i]));
                RUN(json_append(o, ", \"label\": "));
                RUN(json_str(o, s->mapq_map->description[i]));
                RUN(json_append(o, i + 1 < JSON_N_BIN ? "},\n" : "}\n"));
        }
        RUN(json_append(o, "],\n\"axes\": {\n\"length_edges\": ["));
        for(int i = 0; i < COLLECT_N_LEN_BINS; i++){
                if(i){
                        RUN(json_char(o, ','));
                }
                RUN(json_u64(o, collect_length_edge(i)));
        }
        RUN(json_append(o, "],\n\"absolute\": {\"bins\": "));
        RUN(json_u64(o, COLLECT_ABS_BINS));
        RUN(json_append(o, ", \"label\": \"position in read (bp)\"},\n"));
        RUN(json_append(o, "\"relative\": {\"bins\": "));
        RUN(json_u64(o, COLLECT_REL_BINS));
        RUN(json_append(o, ", \"label\": \"position in read (% of length)\"},\n"));
        RUN(json_append(o, "\"identity_bins\": "));
        RUN(json_u64(o, COLLECT_N_IDENT));
        RUN(json_append(o, ",\n\"quality_error_scale\": "));
        RUN(json_u64(o, COLLECT_QUAL_ONE));
        RUN(json_append(o, "\n},\n\"read1\": "));
        RUN(json_v2_read(o, s, 0));
        if(s->n_read2){
                RUN(json_append(o, ",\n\"read2\": "));
                RUN(json_v2_read(o, s, 1));
        }
        RUN(json_append(o, "\n}"));
        return OK;
ERROR:
        return FAIL;
}

#define V2_U(name, val) do{ RUN(json_append(o, "\"" name "\": ")); RUN(json_u64(o, (val))); RUN(json_append(o, ",\n")); }while(0)
#define V2_B(name, val) do{ RUN(json_append(o, "\"" name "\": ")); RUN(json_append(o, (val) ? "true" : "false")); RUN(json_append(o, ",\n")); }while(0)
#define V2_D(name, val) do{ RUN(json_append(o, "\"" name "\": ")); RUN(json_double(o, (val))); RUN(json_append(o, ",\n")); }while(0)

int json_v2_summary(struct jout *o, struct stat_collection *s)
{
        RUN(json_append(o, "\"summary\": {\n"));
        V2_U("total_reads", s->n_total);
        V2_U("read1_reads", s->n_read1);
        V2_U("read2_reads", s->n_read2);
        V2_U("paired", s->n_paired);
        V2_U("properly_paired", s->n_proper_paired);
        V2_U("secondary", s->n_secondary);
        V2_U("supplementary", s->n_supplementary);
        V2_B("is_aligned", s->is_aligned);
        V2_B("is_partial_report", s->is_partial_report);
        V2_U("total_bases", s->total_bases);
        V2_U("min_length", s->n_total ? s->min_length : 0);
        V2_U("max_length", s->max_length);
        V2_U("n50", s->n50);
        V2_U("median_length", s->median_length);
        V2_D("mean_q", s->mean_q);
        V2_D("pct_mapped", s->pct_mapped);
        V2_D("pct_mapq_ge30", s->pct_mapq_ge30);
        V2_D("pct_clipped", s->pct_clipped);
        V2_D("median_identity", s->median_identity);
        RUN(json_append(o, "\"layout\": "));
        RUN(json_str(o, s->layout == COLLECT_LAYOUT_LONG ? "long" : "short"));
        RUN(json_append(o, ",\n"));
        V2_U("reads_without_md", s->reads_without_md);
        V2_U("identity_reads", s->identity_reads);
        V2_U("identity_na", s->identity_na);
        V2_U("bases_beyond_position_axis", s->bases_beyond_position_axis);
        RUN(json_append(o, "\"mapq_counts\": {\"read1\": "));
        RUN(json_bin_counts(o, s->basic_nums[0]));
        if(s->n_read2){
                RUN(json_append(o, ", \"read2\": "));
                RUN(json_bin_counts(o, s->basic_nums[1]));
        }
        RUN(json_append(o, "}\n}"));
        return OK;
ERROR:
        return FAIL;
}

/* {"mapq_ge30": [...], ...}: for class m the n values at base + m * class_stride + i * stride,
 * trailing zeros trimmed. */
int json_v2_per_class(struct jout *o, const char *base, size_t class_stride, size_t stride, int n)
{
        RUN(json_append(o, "{"));
        for(int m = 0; m < JSON_N_BIN; m++){
                RUN(json_key(o, bin_id[m]));
                RUN(json_strided(o, base + (size_t) m * class_stride, stride, n));
                if(m + 1 < JSON_N_BIN){
                        RUN(json_append(o, ", "));
                }
        }
        RUN(json_append(o, "}"));
        return OK;
ERROR:
        return FAIL;
}

/* {"mapq_ge30": {"A": [...], "C": [...], ...}, ...}: as json_v2_per_class, with the five base
 * counters (A C G T N, consecutive uint64_t) of each element at base */
int json_v2_per_class_base(struct jout *o, const char *base, size_t class_stride, size_t stride, int n)
{
        RUN(json_append(o, "{"));
        for(int m = 0; m < JSON_N_BIN; m++){
                RUN(json_key(o, bin_id[m]));
                RUN(json_append(o, "{"));
                for(int k = 0; k < JSON_N_BASE; k++){
                        RUN(json_key(o, base_name[k]));
                        RUN(json_strided(o, base + (size_t) m * class_stride + (size_t) k * sizeof(uint64_t), stride, n));
                        if(k + 1 < JSON_N_BASE){
                                RUN(json_append(o, ", "));
                        }
                }
                RUN(json_append(o, m + 1 < JSON_N_BIN ? "}, " : "}"));
        }
        RUN(json_append(o, "}"));
        return OK;
ERROR:
        return FAIL;
}

/* bins = acc->pos_*[r][0]; the classes follow at a stride of n_bins bins. mismatch_base (added
 * after 2.4.0, format_version 2) follows mismatches. */
int json_v2_pos_block(struct jout *o, const struct collect_pos_bin *bins, int n_bins)
{
        const char *b = (const char *) bins;
        size_t cs = (size_t) n_bins * sizeof(struct collect_pos_bin);
        size_t st = sizeof(struct collect_pos_bin);
        static const char *const names[] = {"qual_bases", "qual_error_sum", "aligned", "aligned_md",
                "mismatches", "insertions", "deletions"};
        static const size_t offs[] = {offsetof(struct collect_pos_bin, qual_bases),
                offsetof(struct collect_pos_bin, qual_error_sum), offsetof(struct collect_pos_bin, aligned),
                offsetof(struct collect_pos_bin, aligned_md), offsetof(struct collect_pos_bin, mismatches),
                offsetof(struct collect_pos_bin, insertions), offsetof(struct collect_pos_bin, deletions)};

        RUN(json_append(o, "{\n"));
        RUN(json_key(o, "composition"));
        RUN(json_v2_per_class_base(o, b + offsetof(struct collect_pos_bin, comp), cs, st, n_bins));
        for(int i = 0; i < 7; i++){
                RUN(json_append(o, ",\n"));
                RUN(json_key(o, names[i]));
                RUN(json_v2_per_class(o, b + offs[i], cs, st, n_bins));
                if(offs[i] == offsetof(struct collect_pos_bin, mismatches)){
                        RUN(json_append(o, ",\n"));
                        RUN(json_key(o, "mismatch_base"));
                        RUN(json_v2_per_class_base(o, b + offsetof(struct collect_pos_bin, mismatch_base), cs, st, n_bins));
                }
        }
        RUN(json_append(o, "\n}"));
        return OK;
ERROR:
        return FAIL;
}

int json_v2_read(struct jout *o, struct stat_collection *s, int r)
{
        struct collect_acc *a = s->acc;
        size_t u = sizeof(uint64_t);

        RUN(json_append(o, "{\n\"n_reads\": "));
        RUN(json_u64(o, r ? s->n_read2 : s->n_read1));
        RUN(json_append(o, ",\n\"length\": {\"reads\": "));
        RUN(json_v2_per_class(o, (const char *) a->len_reads[r], COLLECT_N_LEN_BINS * u, u, COLLECT_N_LEN_BINS));
        RUN(json_append(o, ", \"bases\": "));
        RUN(json_v2_per_class(o, (const char *) a->len_bases[r], COLLECT_N_LEN_BINS * u, u, COLLECT_N_LEN_BINS));
        RUN(json_append(o, "},\n\"read_mean_q\": "));
        RUN(json_v2_per_class(o, (const char *) a->read_q[r], COLLECT_N_Q * u, u, COLLECT_N_Q));
        RUN(json_append(o, ",\n\"identity\": "));
        RUN(json_v2_per_class(o, (const char *) a->identity[r], COLLECT_N_IDENT * u, u, COLLECT_N_IDENT));
        RUN(json_append(o, ",\n\"absolute\": "));
        RUN(json_v2_pos_block(o, a->pos_abs[r][0], COLLECT_ABS_BINS));
        RUN(json_append(o, ",\n\"relative\": "));
        RUN(json_v2_pos_block(o, a->pos_rel[r][0], COLLECT_REL_BINS));
        RUN(json_append(o, "\n}"));
        return OK;
ERROR:
        return FAIL;
}

/* n uint64 values at base + i * stride, trailing zeros trimmed */
int json_strided(struct jout *o, const char *base, size_t stride, int n)
{
        int w = 0;
        uint64_t v;
        for(int i = n - 1; i >= 0; i--){
                memcpy(&v, base + (size_t) i * stride, sizeof(v));
                if(v){
                        w = i + 1;
                        break;
                }
        }
        RUN(json_char(o, '['));
        for(int i = 0; i < w; i++){
                if(i){
                        RUN(json_char(o, ','));
                }
                memcpy(&v, base + (size_t) i * stride, sizeof(v));
                RUN(json_u64(o, v));
        }
        RUN(json_char(o, ']'));
        return OK;
ERROR:
        return FAIL;
}

/* floating point values are printed %.4f, a value that does not apply as null */
int json_double(struct jout *o, double v)
{
        char buf[64];
        if(COLLECT_IS_NA(v)){
                RUN(json_append(o, "null"));
        }else{
                snprintf(buf, 64, "%.4f", v);
                RUN(json_append(o, buf));
        }
        return OK;
ERROR:
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
