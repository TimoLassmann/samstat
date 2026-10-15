/* Unit test for the v2 collector (src/collect/collect.c, standard-plots.md section 3).
 *
 * Files are read the way samstat.c reads them (reader, parse_alignment, fix_orientation,
 * collect_stats) and the accumulators and the finalised summary values are compared with
 * numbers that do not come from samstat: the ground truth of the fixture generator
 * (tests/data/longread/README.md) and values worked out by hand.
 *
 *   merge        a zero-read collection is the identity element; merge is order independent
 *                and independent of how the records are split over collections
 *   tables       E[0] = 2^30, length edges and display edges (every display bin is a union
 *                of storage bins)
 *   bins         relative and absolute position bins, bases beyond the absolute axis
 *   per record   quality is counted per record, never decided from the first chunk
 *   mean Q       128-bit error sum above UINT64_MAX
 *   fixtures     summary values of the long-read fixtures
 */
#include "tldlite.h"
#include "../../src/htsinterface/htsglue.h"
#include "../../src/sambamparse/sam_bam_parse.h"
#include "../../src/collect/collect.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef COLLECT_TEST_DATA_DIR
#error "COLLECT_TEST_DATA_DIR (tests/data) must be defined"
#endif

static int n_checks = 0;
static int n_failed = 0;

#define CHECK(COND, ...) do{                                            \
                n_checks++;                                             \
                if(!(COND)){                                            \
                        n_failed++;                                     \
                        fprintf(stderr, "FAIL (%s:%d): %s: ", __FILE__, __LINE__, #COND); \
                        fprintf(stderr, __VA_ARGS__);                   \
                        fprintf(stderr, "\n");                          \
                }                                                       \
        }while(0)

#define CHECK_U(WHAT, GOT, WANT) do{                                    \
                unsigned long long g_ = (unsigned long long)(GOT);      \
                unsigned long long w_ = (unsigned long long)(WANT);     \
                CHECK(g_ == w_, "%s: %s is %llu, expected %llu", ctx, WHAT, g_, w_); \
        }while(0)

/* decimals are compared to the 4 digits the README gives; NA is expected as "NA" */
#define CHECK_D(WHAT, GOT, WANT) do{                                    \
                double g_ = (GOT);                                      \
                double w_ = (WANT);                                     \
                if(w_ == NCD){                                          \
                        n_checks++;                                     \
                }else if(COLLECT_IS_NA(w_)){                            \
                        CHECK(COLLECT_IS_NA(g_), "%s: %s is %.4f, expected NA", ctx, WHAT, g_); \
                }else{                                                  \
                        CHECK(!COLLECT_IS_NA(g_) && fabs(g_ - w_) < 6e-5, "%s: %s is %.4f, expected %.4f", ctx, WHAT, g_, w_); \
                }                                                       \
        }while(0)

#define NA COLLECT_NA
#define NCD (-1.0)      /* decimal not checked */

/* --- reading ---------------------------------------------------------------------- */

/* Read <path> into <n> collections (chunk k into collection k % n). The collections are
 * allocated here. has_sq is passed on from the reader. */
static int read_into(const char *path, struct stat_collection **cols, int n, int batch_records)
{
        struct sam_bam_file *f = NULL;
        struct tl_seq_buffer *sb = NULL;
        struct reader_opts o;
        int chunk = 0;

        o.pool_threads = 0;
        o.batch_records = (uint64_t) batch_records;
        o.batch_bases = 8000000;
        o.max_primary = UINT64_MAX;
        o.max_read_len = 0;

        for(int i = 0; i < n; i++){
                cols[i] = NULL;
                RUN(stat_collection_alloc(&cols[i]));
        }
        RUN(open_reads(&f, (char *) path, &o));
        RUN(alloc_tl_seq_buffer(&sb, batch_records));
        if(!f->is_fastx){
                RUN(add_aln_data(sb));
        }

        while(1){
                RUN(read_bam_chunk(f, sb));
                sb->L = TL_SEQ_BUFFER_DNA;
                sb->base_quality_offset = 33;
                if(sb->num_seq == 0){
                        break;
                }
                if(!f->is_fastx){
                        for(int i = 0; i < sb->num_seq; i++){
                                RUN(parse_alignment(sb->sequences[i]));
                                RUN(fix_orientation(sb->sequences[i]));
                        }
                }
                RUN(collect_stats(sb, cols[chunk % n]));
                chunk++;
                if(!f->is_fastx){
                        clear_aln_data(sb);
                }
                reset_tl_seq_buffer(sb);
        }
        for(int i = 0; i < n; i++){
                cols[i]->has_sq = f->has_sq;
        }
        cols[0]->n_secondary = f->n_secondary;
        cols[0]->n_supplementary = f->n_supplementary;
        RUN(close_bam(f));
        f = NULL;
        free_aln_data(sb);
        free_tl_seq_buffer(sb);
        return OK;
ERROR:
        return FAIL;
}

static struct stat_collection *read_one(const char *path)
{
        struct stat_collection *cols[1];
        if(read_into(path, cols, 1, 4096) != OK){
                fprintf(stderr, "FAIL: cannot read %s\n", path);
                n_failed++;
                return NULL;
        }
        if(stat_collection_finalise(cols[0]) != OK){
                fprintf(stderr, "FAIL: finalise failed for %s\n", path);
                n_failed++;
                stat_collection_free(cols[0]);
                return NULL;
        }
        return cols[0];
}

static void fixture_path(char *buf, size_t n, const char *name)
{
        snprintf(buf, n, "%s/%s", COLLECT_TEST_DATA_DIR, name);
}

/* --- comparing collections ----------------------------------------------------------- */

static int same_double(double a, double b)
{
        if(COLLECT_IS_NA(a) || COLLECT_IS_NA(b)){
                return COLLECT_IS_NA(a) && COLLECT_IS_NA(b);
        }
        return a == b;
}

/* both finalised */
static int same_collection(const struct stat_collection *a, const struct stat_collection *b, const char *ctx)
{
        int ok = 1;
#define SAME(F) if(a->F != b->F){ fprintf(stderr, "  %s: %s differs: %llu vs %llu\n", ctx, #F, (unsigned long long) a->F, (unsigned long long) b->F); ok = 0; }
        SAME(n_total) SAME(n_read1) SAME(n_read2) SAME(n_paired) SAME(n_proper_paired)
        SAME(n_secondary) SAME(n_supplementary) SAME(total_bases) SAME(min_length) SAME(max_length)
        SAME(mapped_bases) SAME(clipped_bases) SAME(qual_bases_all) SAME(reads_without_md)
        SAME(identity_reads) SAME(identity_na) SAME(bases_beyond_position_axis) SAME(has_sq)
        SAME(is_partial_report) SAME(is_aligned) SAME(layout) SAME(n50) SAME(median_length)
#undef SAME
        if(a->qual_error_sum_all != b->qual_error_sum_all){
                fprintf(stderr, "  %s: qual_error_sum_all differs\n", ctx);
                ok = 0;
        }
        if(!same_double(a->mean_q, b->mean_q) || !same_double(a->pct_mapped, b->pct_mapped)
           || !same_double(a->pct_mapq_ge30, b->pct_mapq_ge30) || !same_double(a->pct_clipped, b->pct_clipped)
           || !same_double(a->median_identity, b->median_identity)){
                fprintf(stderr, "  %s: summary percentages differ\n", ctx);
                ok = 0;
        }
        for(int r = 0; r < COLLECT_N_MATE; r++){
                for(int m = 0; m < COLLECT_N_MAPQ; m++){
                        if(a->basic_nums[r][m] != b->basic_nums[r][m]){
                                fprintf(stderr, "  %s: basic_nums[%d][%d] differs\n", ctx, r, m);
                                ok = 0;
                        }
                }
        }
        if(memcmp(a->acc, b->acc, sizeof(struct collect_acc)) != 0){
                fprintf(stderr, "  %s: accumulators differ\n", ctx);
                ok = 0;
        }
        return ok;
}

/* --- merge ---------------------------------------------------------------------------- */

static int test_merge(void)
{
        const char *files[] = {"longread/ont_like.sam", "longread/single_len.sam", "paired_small.sam", "longread/long.fq"};
        char path[1024];

        for(size_t k = 0; k < sizeof(files) / sizeof(files[0]); k++){
                struct stat_collection *whole = NULL;
                struct stat_collection *filled = NULL;
                struct stat_collection *zero = NULL;
                struct stat_collection *zero2 = NULL;
                struct stat_collection *parts[3];
                struct stat_collection *parts2[3];
                struct stat_collection *p3[1];
                char ctx[1100];

                fixture_path(path, sizeof(path), files[k]);
                snprintf(ctx, sizeof(ctx), "merge %s", files[k]);

                whole = read_one(path);
                filled = read_one(path);
                CHECK(whole != NULL && filled != NULL, "%s: cannot read", ctx);
                if(!whole || !filled){
                        continue;
                }
                RUN(stat_collection_alloc(&zero));
                RUN(stat_collection_alloc(&zero2));

                /* filled + zero == filled; zero + filled == filled */
                RUN(stat_collection_merge(filled, zero));
                RUN(stat_collection_finalise(filled));
                CHECK(same_collection(filled, whole, ctx), "%s: filled + zero differs from filled", ctx);

                RUN(stat_collection_merge(zero2, whole));
                RUN(stat_collection_finalise(zero2));
                CHECK(same_collection(zero2, whole, ctx), "%s: zero + filled differs from filled", ctx);

                /* records split over three collections in chunks of 3: merge in two orders */
                RUN(read_into(path, parts, 3, 3));
                RUN(read_into(path, parts2, 3, 3));
                RUN(stat_collection_merge(parts[0], parts[1]));
                RUN(stat_collection_merge(parts[0], parts[2]));
                RUN(stat_collection_merge(parts2[2], parts2[1]));
                RUN(stat_collection_merge(parts2[2], parts2[0]));
                RUN(stat_collection_finalise(parts[0]));
                RUN(stat_collection_finalise(parts2[2]));
                /* the reader counters are set on one collection only: both orders carry them */
                CHECK(same_collection(parts[0], parts2[2], ctx), "%s: merge order matters", ctx);
                CHECK(same_collection(parts[0], whole, ctx), "%s: merged parts differ from the whole file", ctx);

                /* one collection, chunks of 3 records: same as chunks of 4096 */
                RUN(read_into(path, p3, 1, 3));
                RUN(stat_collection_finalise(p3[0]));
                CHECK(same_collection(p3[0], whole, ctx), "%s: the batch size changes the result", ctx);

                for(int i = 0; i < 3; i++){
                        stat_collection_free(parts[i]);
                        stat_collection_free(parts2[i]);
                }
                stat_collection_free(p3[0]);
                stat_collection_free(zero);
                stat_collection_free(zero2);
                stat_collection_free(whole);
                stat_collection_free(filled);
        }

        /* a zero-read collection on its own */
        {
                struct stat_collection *z = NULL;
                const char *ctx = "empty";
                RUN(stat_collection_alloc(&z));
                CHECK_U("min_length before finalise", z->min_length, UINT64_MAX);
                RUN(stat_collection_finalise(z));
                CHECK_U("total", z->n_total, 0);
                CHECK_U("min_length", z->min_length, 0);
                CHECK_U("n50", z->n50, 0);
                CHECK_U("median_length", z->median_length, 0);
                CHECK_U("is_aligned", z->is_aligned, 0);
                CHECK_D("mean_q", z->mean_q, NA);
                CHECK_D("pct_mapped", z->pct_mapped, NA);
                CHECK_D("pct_clipped", z->pct_clipped, NA);
                CHECK_D("median_identity", z->median_identity, NA);
                stat_collection_free(z);
        }
        return OK;
ERROR:
        n_failed++;
        return FAIL;
}

/* --- tables --------------------------------------------------------------------------- */

static int test_tables(void)
{
        const char *ctx = "tables";
        int n_display = collect_n_display_bins();

        CHECK_U("E[0]", collect_qual_error(0), 1073741824ull);
        CHECK_U("E[0] is 2^30", collect_qual_error(0), 1ull << 30);
        CHECK_U("E[10] = round(0.1 * 2^30)", collect_qual_error(10), 107374182ull);
        CHECK_U("E[20] = round(0.01 * 2^30)", collect_qual_error(20), 10737418ull);
        CHECK_U("E[93]", collect_qual_error(93), 1);
        for(int q = 1; q < COLLECT_N_Q; q++){
                CHECK(collect_qual_error(q) <= collect_qual_error(q - 1), "E[%d] grows", q);
                CHECK(collect_qual_error(q) >= 1, "E[%d] is 0", q);
        }

        /* storage bins: 0..999 exact, then the 1% grid */
        CHECK_U("storage bins", COLLECT_N_LEN_BINS, 1922);
        for(int i = 0; i < COLLECT_LEN_EXACT; i++){
                CHECK(collect_length_edge(i) == (uint64_t) i && collect_length_bin((uint64_t) i) == i, "exact length %d", i);
        }
        CHECK_U("edge 1000", collect_length_edge(1000), 1000);
        CHECK_U("edge 1001", collect_length_edge(1001), 1010);
        CHECK_U("edge 1002", collect_length_edge(1002), 1021);
        CHECK_U("last edge", collect_length_edge(COLLECT_N_LEN_BINS - 1), 10015438);
        for(int b = 1; b < COLLECT_N_LEN_BINS; b++){
                uint64_t e = collect_length_edge(b);
                uint64_t prev = collect_length_edge(b - 1);
                /* each bin is at most 1% wide (plus rounding) above 1 kb, and non-empty */
                CHECK(e > prev, "edge %d not increasing", b);
                CHECK(collect_length_bin(e) == b, "edge %llu is in bin %d", (unsigned long long) e, collect_length_bin(e));
                CHECK(collect_length_bin(e - 1) == b - 1, "edge %llu - 1 is in bin %d", (unsigned long long) e, collect_length_bin(e - 1));
        }
        CHECK_U("length 1 Gb is in the overflow bin", collect_length_bin(1000000000ull), COLLECT_N_LEN_BINS - 1);

        /* display bins: 0, g(1) below 1000, then the storage edges from 1000 */
        CHECK_U("display bins", n_display, 1217);
        CHECK_U("display edge 0", collect_display_edge(0), 0);
        CHECK_U("display edge 1", collect_display_edge(1), 1);
        CHECK_U("display edge 100", collect_display_edge(100), 100);
        CHECK_U("display edge 101", collect_display_edge(101), 101);
        CHECK_U("display edge 102", collect_display_edge(102), 103);
        CHECK_U("last display edge", collect_display_edge(n_display - 1), 10015438);
        for(int d = 1; d < n_display; d++){
                CHECK(collect_display_edge(d) > collect_display_edge(d - 1), "display edge %d not increasing", d);
        }
        /* every storage edge from 1000 on is a display edge ... */
        {
                int d = 0;
                for(int b = COLLECT_LEN_EXACT; b < COLLECT_N_LEN_BINS; b++){
                        while(d < n_display && collect_display_edge(d) < collect_length_edge(b)){
                                d++;
                        }
                        CHECK(d < n_display && collect_display_edge(d) == collect_length_edge(b),
                              "storage edge %llu is not a display edge", (unsigned long long) collect_length_edge(b));
                }
        }
        /* ... and below 1000 every display edge is a storage edge, so each display bin is a
         * union of whole storage bins: the display bin never changes inside a display bin */
        for(int d = 0; d < n_display; d++){
                uint64_t e = collect_display_edge(d);
                CHECK(e < 1000 || collect_length_edge(collect_length_bin(e)) == e, "display edge %llu is not a storage edge", (unsigned long long) e);
        }
        {
                int last = 0;
                for(int b = 0; b < COLLECT_N_LEN_BINS; b++){
                        int d = collect_display_bin(b);
                        CHECK(d >= last && d - last <= 1, "display bin of storage bin %d jumps from %d to %d", b, last, d);
                        CHECK(collect_display_edge(d) <= collect_length_edge(b), "display bin %d starts above storage bin %d", d, b);
                        if(d + 1 < n_display){
                                CHECK(collect_display_edge(d + 1) > collect_length_edge(b), "display bin %d ends below storage bin %d", d, b);
                        }
                        last = d;
                }
        }
        return OK;
}

/* --- position bins ---------------------------------------------------------------------- */

#define SAM_HEADER "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n"

/* One unmapped read of <len> bases A, Q40. Returns the collection (finalised). */
static struct stat_collection *read_synthetic(const char *name, const char *records, int batch_records, int *ok)
{
        struct stat_collection *cols[1];
        FILE *fp = fopen(name, "w");

        *ok = 0;
        if(!fp){
                fprintf(stderr, "FAIL: cannot write %s\n", name);
                n_failed++;
                return NULL;
        }
        fputs(SAM_HEADER, fp);
        fputs(records, fp);
        if(fclose(fp) != 0){
                n_failed++;
                return NULL;
        }
        if(read_into(name, cols, 1, batch_records) != OK || stat_collection_finalise(cols[0]) != OK){
                fprintf(stderr, "FAIL: cannot read %s\n", name);
                n_failed++;
                remove(name);
                return NULL;
        }
        remove(name);
        *ok = 1;
        return cols[0];
}

static int test_position_bins(void)
{
        const int lens[] = {1, 7, 37, 99, 100, 101, 999, 1000, 1001, 2500};
        char name[64] = "unit_collect_pos.sam";
        char ctx[64];

        for(size_t k = 0; k < sizeof(lens) / sizeof(lens[0]); k++){
                const int L = lens[k];
                char *rec = malloc((size_t) 2 * L + 128);
                struct stat_collection *s = NULL;
                int ok = 0;
                uint64_t expect_rel[COLLECT_REL_BINS];

                snprintf(ctx, sizeof(ctx), "position bins L=%d", L);
                if(!rec){
                        n_failed++;
                        continue;
                }
                /* an unmapped read: the whole read is the 'aligned part' */
                {
                        int n = snprintf(rec, 64, "r1\t4\t*\t0\t0\t*\t*\t0\t0\t");
                        memset(rec + n, 'A', (size_t) L);
                        n += L;
                        rec[n++] = '\t';
                        memset(rec + n, 'I', (size_t) L);
                        n += L;
                        rec[n++] = '\n';
                        rec[n] = 0;
                }
                s = read_synthetic(name, rec, 4096, &ok);
                free(rec);
                if(!ok){
                        continue;
                }
                for(int b = 0; b < COLLECT_REL_BINS; b++){
                        expect_rel[b] = 0;
                }
                for(int i = 0; i < L; i++){
                        expect_rel[(100 * i) / L]++;
                }
                for(int b = 0; b < COLLECT_REL_BINS; b++){
                        const struct collect_pos_bin *pb = &s->acc->pos_rel[0][COLLECT_MAPQ_UNMAPPED][b];
                        CHECK(pb->comp[0] == expect_rel[b], "%s: relative bin %d has %llu A, expected %llu", ctx, b, (unsigned long long) pb->comp[0], (unsigned long long) expect_rel[b]);
                        CHECK(pb->qual_bases == expect_rel[b], "%s: relative bin %d has %llu qualities", ctx, b, (unsigned long long) pb->qual_bases);
                        CHECK(pb->qual_error_sum == expect_rel[b] * collect_qual_error(40), "%s: relative bin %d error sum", ctx, b);
                }
                for(int i = 0; i < COLLECT_ABS_BINS; i++){
                        const struct collect_pos_bin *pb = &s->acc->pos_abs[0][COLLECT_MAPQ_UNMAPPED][i];
                        uint64_t want = (i < L) ? 1 : 0;
                        CHECK(pb->comp[0] == want && pb->qual_bases == want, "%s: absolute bin %d has %llu A, expected %llu", ctx, i, (unsigned long long) pb->comp[0], (unsigned long long) want);
                }
                CHECK_U("bases beyond the absolute axis", s->bases_beyond_position_axis, (L > 1000) ? L - 1000 : 0);
                CHECK_U("total_bases", s->total_bases, L);
                CHECK_U("qual_bases_all", s->qual_bases_all, L);
                CHECK_D("mean_q", s->mean_q, 40.0);
                /* never aligned: no aligned bases anywhere */
                {
                        uint64_t al = 0;
                        for(int b = 0; b < COLLECT_REL_BINS; b++){
                                al += s->acc->pos_rel[0][COLLECT_MAPQ_UNMAPPED][b].aligned;
                        }
                        CHECK_U("aligned bases of an unmapped read", al, 0);
                }
                stat_collection_free(s);
        }
        return OK;
}

/* --- per-record decisions ----------------------------------------------------------------- */

static int test_per_record(void)
{
        const char *ctx = "per record";
        char name[64] = "unit_collect_mixed.sam";
        /* record 1 has no quality, 2 and 3 have Q40 (I): with one record per chunk the old
         * first-chunk decision would have dropped all qualities */
        const char *recs =
                "r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n"
                "r2\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
                "r3\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";
        int ok = 0;
        struct stat_collection *s = read_synthetic(name, recs, 1, &ok);
        uint64_t n_q = 0;

        if(!ok){
                return OK;
        }
        CHECK_U("total_reads", s->n_total, 3);
        CHECK_U("qual_bases_all", s->qual_bases_all, 20);
        CHECK_D("mean_q", s->mean_q, 40.0);
        for(int q = 0; q < COLLECT_N_Q; q++){
                n_q += s->acc->read_q[0][COLLECT_MAPQ_GE30][q];
        }
        CHECK_U("reads with a mean Q", n_q, 2);
        CHECK_U("reads in Q bin 39 or 40", s->acc->read_q[0][COLLECT_MAPQ_GE30][39] + s->acc->read_q[0][COLLECT_MAPQ_GE30][40], 2);
        /* the position bins hold the quality of the two reads and the bases of all three */
        CHECK_U("bases at position 0", s->acc->pos_abs[0][COLLECT_MAPQ_GE30][0].comp[0] + s->acc->pos_abs[0][COLLECT_MAPQ_GE30][0].comp[1]
                + s->acc->pos_abs[0][COLLECT_MAPQ_GE30][0].comp[2] + s->acc->pos_abs[0][COLLECT_MAPQ_GE30][0].comp[3], 3);
        CHECK_U("qualities at position 0", s->acc->pos_abs[0][COLLECT_MAPQ_GE30][0].qual_bases, 2);
        CHECK_U("reads without MD", s->reads_without_md, 3);
        stat_collection_free(s);
        return OK;
}

/* --- identity of a hand-computed CIGAR ------------------------------------------------------ */

static int test_identity_bins(void)
{
        const char *ctx = "identity_indel.sam";
        char path[1024];
        struct stat_collection *s;
        uint64_t sum = 0;

        fixture_path(path, sizeof(path), "longread/identity_indel.sam");
        s = read_one(path);
        if(!s){
                return OK;
        }
        /* id_ins1 99/100 = 99.0% (bin 990), id_del2_mm1 99/101 = 98.0198% (bin 980),
         * id_perfect 100% (bin 1000), id_ins3 97/98 = 98.9796% (bin 989) */
        {
                int want_bins[] = {980, 989, 990, 1000};
                for(int b = 0; b < COLLECT_N_IDENT; b++){
                        for(int m = 0; m < COLLECT_N_MAPQ; m++){
                                sum += s->acc->identity[0][m][b];
                        }
                }
                CHECK_U("binned reads", sum, 4);
                for(int k = 0; k < 4; k++){
                        uint64_t n = 0;
                        for(int m = 0; m < COLLECT_N_MAPQ; m++){
                                n += s->acc->identity[0][m][want_bins[k]];
                        }
                        CHECK_U("reads in identity bin", n, 1);
                }
        }
        CHECK_D("median identity (lower median: bin 989)", s->median_identity, 98.9);
        stat_collection_free(s);
        return OK;
}

/* --- mean Q above UINT64_MAX -------------------------------------------------------------------- */

static int test_mean_q_128(void)
{
        struct stat_collection *s = NULL;
        const uint64_t bases = 1ull << 47;      /* 1.4e14 bases */
        const uint64_t e10 = collect_qual_error(10);    /* Q10: 0.1 */

        if(stat_collection_alloc(&s) != OK){
                n_failed++;
                return FAIL;
        }
        s->qual_bases_all = bases;
        s->qual_error_sum_all = (collect_u128) e10 * bases;
        CHECK(s->qual_error_sum_all > (collect_u128) UINT64_MAX, "the synthetic sum does not exceed UINT64_MAX");
        CHECK(stat_collection_finalise(s) == OK, "finalise");
        CHECK(!COLLECT_IS_NA(s->mean_q) && fabs(s->mean_q - 10.0) < 1e-6, "mean_q is %.9f, expected 10", s->mean_q);

        /* Q0 everywhere: sum is 2^30 * 2^47 = 2^77 */
        s->qual_error_sum_all = (collect_u128) collect_qual_error(0) * bases;
        CHECK(stat_collection_finalise(s) == OK, "finalise");
        CHECK(fabs(s->mean_q) < 1e-9, "mean_q is %.9f, expected 0", s->mean_q);
        CHECK(!signbit(s->mean_q), "mean_q of Q0 is -0.0, it prints as \"-0.0\"");
        CHECK(!signbit(collect_mean_q((collect_u128) collect_qual_error(0) * 1000, 1000)), "collect_mean_q of Q0 is -0.0");

        CHECK(fabs(collect_mean_q((collect_u128) collect_qual_error(30) * 1000, 1000) - 30.0) < 1e-6, "collect_mean_q Q30");
        CHECK(COLLECT_IS_NA(collect_mean_q(0, 0)), "collect_mean_q of no bases is NA");
        stat_collection_free(s);
        return OK;
}

/* --- fixtures: values from tests/data/longread/README.md ------------------------------------- */

struct fixture_expect {
        const char *file;
        uint64_t total, r1, r2, paired, proper, secondary, supplementary;
        uint64_t bases, min_len, max_len, n50, median;
        int layout;                   /* COLLECT_LAYOUT_* */
        int is_aligned;
        uint64_t mapq[4];             /* >= 30, 10-29, < 10, unmapped (both mates) */
        double mean_q, pct_mapped, pct_ge30, pct_clipped, median_identity;
        int64_t identity_reads, identity_na, without_md;   /* -1: not checked */
        int64_t mismatch, ins, del, aligned, aligned_md, clipped;   /* sums over the relative axis, -1: not checked */
};

#define NC (-1)

static const struct fixture_expect fixtures[] = {
        {"longread/ont_like.sam", 24, 24, 0, 0, 0, 0, 0, 36990, 538, 3442, 2076, 1283, COLLECT_LAYOUT_LONG, 1,
         {12, 4, 6, 2}, 10.4577, 91.6667, 50.0, 5.5942, 92.9, 22, 0, 0, 949, 533, 716, 31338, 31338, 1857},
        {"longread/ont_like.bam", 25, 25, 0, 0, 0, 0, 0, 89490, 538, 52500, 52310, 1422, COLLECT_LAYOUT_LONG, 1,
         {13, 4, 6, 2}, 13.6586, 92.0, 52.0, 2.1670, 92.9, 23, 0, 0, 949, 18033, 18216, 83838, 83838, 1857},
        {"longread/supp_chain.sam", 2, 2, 0, 0, 0, 1, 2, 3000, 1000, 2000, 1993, 1000, COLLECT_LAYOUT_LONG, 1,
         {2, 0, 0, 0}, 15.1028, 100.0, 100.0, 33.3333, 99.9, 2, 0, NC, 1, 0, 0, 2000, 2000, 1000},
        {"longread/hifi_like.sam", 4, 4, 0, 0, 0, 0, 0, 18000, 3000, 6000, 4757, 4173, COLLECT_LAYOUT_LONG, 1,
         {3, 1, 0, 0}, 37.0198, 100.0, 75.0, 0.0, 99.7, 4, 0, 0, 22, 14, 12, 18000, 18000, 0},
        {"longread/long.fq", 12, 12, 0, 0, 0, 0, 0, 19391, 1, 5000, 2595, 1193, COLLECT_LAYOUT_LONG, 0,
         {0, 0, 0, 12}, 11.8917, NAN, NAN, NAN, NAN, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {"longread/eqx_nomd.sam", 5, 5, 0, 0, 0, 0, 0, 2530, 300, 700, 600, 510, COLLECT_LAYOUT_SHORT, 1,
         {4, 1, 0, 0}, 13.8782, 100.0, 80.0, 2.6087, 95.0, 5, 0, 0, 81, 22, 19, 2464, 2464, 66},
        {"longread/nonm.sam", 5, 5, 0, 0, 0, 0, 0, 2300, 300, 650, 520, 450, COLLECT_LAYOUT_SHORT, 1,
         {4, 0, 1, 0}, 15.0286, 100.0, 80.0, 0.0, NAN, 0, 5, 0, 32, 17, 26, 2300, 2300, 0},
        {"longread/ubam_like.sam", 6, 6, 0, 0, 0, 0, 0, 4100, 300, 1100, 750, 600, COLLECT_LAYOUT_SHORT, 0,
         {0, 0, 0, 6}, 12.1286, NAN, NAN, NAN, NAN, NC, NC, NC, NC, NC, NC, NC, NC, NC},
        {"longread/sq_nomap.sam", 6, 6, 0, 0, 0, 0, 0, 3820, 300, 1000, 700, 550, COLLECT_LAYOUT_SHORT, 1,
         {0, 0, 0, 6}, 12.0356, 0.0, 0.0, NAN, NAN, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {"longread/mapq255.sam", 7, 7, 0, 0, 0, 0, 0, 3850, 320, 760, 640, 560, COLLECT_LAYOUT_SHORT, 1,
         {6, 0, 1, 0}, 28.0353, 100.0, 85.7143, 0.0, 96.7, 7, 0, 0, 80, 20, 22, 3850, 3850, 0},
        {"longread/single_len.sam", 20, 10, 10, 20, 20, 0, 0, 3000, 150, 150, 150, 150, COLLECT_LAYOUT_SHORT, 1,
         {14, 4, 2, 0}, 30.0519, 100.0, 70.0, 0.0, 99.3, 20, 0, 0, 24, 0, 1, 3000, 3000, 0},
        {"longread/n50_five.sam", 5, NC, NC, NC, NC, NC, NC, 1600, NC, NC, 400, 300, COLLECT_LAYOUT_SHORT, 1,
         {5, NC, NC, NC}, NCD, NCD, NCD, NCD, 100.0, NC, NC, NC, 0, NC, NC, NC, NC, NC},
        {"longread/identity_indel.sam", 4, NC, NC, NC, NC, NC, NC, 400, NC, NC, NC, NC, COLLECT_LAYOUT_SHORT, 1,
         {NC, NC, NC, NC}, NCD, NCD, NCD, NCD, 98.9, 4, 0, NC, 1, 2, 1, 400, NC, NC},
        {"longread/softclip_pct.sam", 5, NC, NC, NC, NC, NC, NC, 500, NC, NC, NC, NC, COLLECT_LAYOUT_SHORT, 1,
         {4, NC, NC, 1}, NCD, 80.0, 80.0, 15.0, NCD, NC, NC, NC, NC, NC, NC, NC, NC, 60},
};

#define MAYBE_U(WHAT, GOT, WANT) if((int64_t)(WANT) != NC){ CHECK_U(WHAT, GOT, WANT); }
/* In the table NAN means "expected not available" for a decimal, NCD "not checked". */

static int test_fixtures(void)
{
        for(size_t k = 0; k < sizeof(fixtures) / sizeof(fixtures[0]); k++){
                const struct fixture_expect *e = &fixtures[k];
                const char *ctx = e->file;
                char path[1024];
                struct stat_collection *s;
                uint64_t mapq[4] = {0, 0, 0, 0};
                int64_t sum_mis = 0, sum_ins = 0, sum_del = 0, sum_al = 0, sum_md = 0;

                fixture_path(path, sizeof(path), e->file);
                s = read_one(path);
                if(!s){
                        continue;
                }
                CHECK_U("total_reads", s->n_total, e->total);
                MAYBE_U("read1_reads", s->n_read1, e->r1);
                MAYBE_U("read2_reads", s->n_read2, e->r2);
                MAYBE_U("paired", s->n_paired, e->paired);
                MAYBE_U("properly_paired", s->n_proper_paired, e->proper);
                MAYBE_U("secondary", s->n_secondary, e->secondary);
                MAYBE_U("supplementary", s->n_supplementary, e->supplementary);
                CHECK_U("total_bases", s->total_bases, e->bases);
                MAYBE_U("min_length", s->min_length, e->min_len);
                MAYBE_U("max_length", s->max_length, e->max_len);
                MAYBE_U("n50", s->n50, e->n50);
                MAYBE_U("median_length", s->median_length, e->median);
                CHECK_U("layout", s->layout, e->layout);
                CHECK_U("is_aligned", s->is_aligned, e->is_aligned);
                for(int r = 0; r < 2; r++){
                        for(int m = 0; m < 4; m++){
                                mapq[m] += s->basic_nums[r][m];
                        }
                }
                for(int m = 0; m < 4; m++){
                        MAYBE_U("MAPQ class count", mapq[m], e->mapq[m]);
                }
                CHECK_D("mean_q", s->mean_q, e->mean_q);
                CHECK_D("pct_mapped", s->pct_mapped, e->pct_mapped);
                CHECK_D("pct_mapq_ge30", s->pct_mapq_ge30, e->pct_ge30);
                CHECK_D("pct_clipped", s->pct_clipped, e->pct_clipped);
                CHECK_D("median_identity", s->median_identity, e->median_identity);
                MAYBE_U("identity_reads", s->identity_reads, e->identity_reads);
                MAYBE_U("identity_na", s->identity_na, e->identity_na);
                MAYBE_U("reads_without_md", s->reads_without_md, e->without_md);
                for(int r = 0; r < COLLECT_N_MATE; r++){
                        for(int m = 0; m < COLLECT_N_MAPQ; m++){
                                for(int b = 0; b < COLLECT_REL_BINS; b++){
                                        const struct collect_pos_bin *pb = &s->acc->pos_rel[r][m][b];
                                        sum_mis += (int64_t) pb->mismatches;
                                        sum_ins += (int64_t) pb->insertions;
                                        sum_del += (int64_t) pb->deletions;
                                        sum_al += (int64_t) pb->aligned;
                                        sum_md += (int64_t) pb->aligned_md;
                                }
                        }
                }
                MAYBE_U("mismatch_bases", sum_mis, e->mismatch);
                MAYBE_U("insertion_events", sum_ins, e->ins);
                MAYBE_U("deletion_events", sum_del, e->del);
                MAYBE_U("aligned_bases", sum_al, e->aligned);
                MAYBE_U("aligned_md_bases", sum_md, e->aligned_md);
                MAYBE_U("clipped_bases", s->clipped_bases, e->clipped);
                /* the relative and the absolute axis hold the same counts while no read is longer than 1000 */
                if(s->max_length <= COLLECT_ABS_BINS){
                        int64_t abs_mis = 0, abs_ins = 0, abs_del = 0;
                        for(int r = 0; r < COLLECT_N_MATE; r++){
                                for(int m = 0; m < COLLECT_N_MAPQ; m++){
                                        for(int b = 0; b < COLLECT_ABS_BINS; b++){
                                                const struct collect_pos_bin *pb = &s->acc->pos_abs[r][m][b];
                                                abs_mis += (int64_t) pb->mismatches;
                                                abs_ins += (int64_t) pb->insertions;
                                                abs_del += (int64_t) pb->deletions;
                                        }
                                }
                        }
                        CHECK(abs_mis == sum_mis && abs_ins == sum_ins && abs_del == sum_del,
                              "%s: absolute axis (%lld %lld %lld) and relative axis (%lld %lld %lld) disagree", ctx,
                              (long long) abs_mis, (long long) abs_ins, (long long) abs_del,
                              (long long) sum_mis, (long long) sum_ins, (long long) sum_del);
                }
                /* summary values are computed from integers: finalise twice gives the same */
                {
                        double mq = s->mean_q;
                        uint64_t n50 = s->n50;
                        CHECK(stat_collection_finalise(s) == OK, "second finalise");
                        CHECK(same_double(mq, s->mean_q) && n50 == s->n50, "%s: finalise is not repeatable", ctx);
                }
                stat_collection_free(s);
        }
        return OK;
}

int main(void)
{
        test_tables();
        test_merge();
        test_position_bins();
        test_per_record();
        test_identity_bins();
        test_mean_q_128();
        test_fixtures();
        if(n_failed){
                fprintf(stderr, "unit_collect: %d of %d checks failed\n", n_failed, n_checks);
                return EXIT_FAILURE;
        }
        fprintf(stderr, "unit_collect: %d checks passed\n", n_checks);
        return EXIT_SUCCESS;
}
