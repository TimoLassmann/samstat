#include "tld.h"

#include "../plot/plot.h"
#include "sam.h"
#include "../htsinterface/htsglue.h"
#include "../param/param.h"

#include "convert_tables.h"

#define COLLECT_IMPORT
#include "collect.h"

#include <math.h>
#include <string.h>

#define MAPQUALBIN_gt30 COLLECT_MAPQ_GE30
#define MAPQUALBIN_lt30 COLLECT_MAPQ_10_29
#define MAPQUALBIN_lt10 COLLECT_MAPQ_LT10
#define MAPQUALBIN_UNMAP COLLECT_MAPQ_UNMAPPED

/* The generated tables and the sizes in collect.h must agree. */
_Static_assert(CT_N_QUAL == COLLECT_N_Q, "quality table size");
_Static_assert(CT_N_LENGTH_EDGES == COLLECT_N_LEN_BINS, "length edge table size");
_Static_assert(sizeof(collect_u128) == 16, "128-bit integers are required");

static int get_mapqual_bins(struct mapqual_bins **map);
static void free_mapqual_bins(struct mapqual_bins *m);

static int collect_stats_end_seqqual(struct tl_seq_buffer *sb, struct stat_collection *s);

static int collect_stats_norm_seqqual(struct tl_seq_buffer *sb, struct stat_collection *s);

#define QUAL_MAX_BIN 41
#define QUAL_MISSING 0xFF

/* Legacy (2.3) collectors only: Phred value of base 'pos' of 'itm' as a plot bin in [0, 41].
   Returns -1 if there is no quality for this base (QUAL '*' is stored as
   0xFF, or the quality string is shorter than the sequence): such bases must
   not contribute to any quality plot. */
static inline int get_base_qual(const struct tl_seq *itm, int pos, int offset)
{
        int q;
        if(pos < 0 || pos >= (int) itm->qual->len){
                return -1;
        }
        if((uint8_t) itm->qual->str[pos] == QUAL_MISSING){
                return -1;
        }
        q = (int)((uint8_t) itm->qual->str[pos]) - offset;
        if(q < 0){
                q = 0;
        }
        if(q > QUAL_MAX_BIN){
                q = QUAL_MAX_BIN;
        }
        return q;
}

/* Clamp a histogram index into [0, len-1]. */
static inline int clamp_idx(int idx, int len)
{
        if(idx >= len){
                idx = len - 1;
        }
        if(idx < 0){
                idx = 0;
        }
        return idx;
}

/* ------------------------------------------------------------------------------------
 * Tables (generated at build time, see convert_tables.c)
 * ------------------------------------------------------------------------------------ */

uint32_t collect_qual_error(int q)
{
        if(q < 0){
                q = 0;
        }
        if(q >= COLLECT_N_Q){
                q = COLLECT_N_Q - 1;
        }
        return qual_error[q];
}

/* Storage bin of a read length: exact below 1000, then the largest edge <= len. Lengths
 * from the last edge on go into the overflow bin. Integer arithmetic only. */
int collect_length_bin(uint64_t len)
{
        int lo = COLLECT_LEN_EXACT;
        int hi = COLLECT_N_LEN_BINS - 1;

        if(len < COLLECT_LEN_EXACT){
                return (int) len;
        }
        if(len >= length_edges[hi]){
                return hi;
        }
        while(lo < hi){
                int mid = (lo + hi + 1) / 2;
                if(length_edges[mid] <= len){
                        lo = mid;
                }else{
                        hi = mid - 1;
                }
        }
        return lo;
}

uint64_t collect_length_edge(int bin)
{
        if(bin < 0){
                bin = 0;
        }
        if(bin >= COLLECT_N_LEN_BINS){
                bin = COLLECT_N_LEN_BINS - 1;
        }
        return length_edges[bin];
}

int collect_n_display_bins(void)
{
        return CT_N_DISPLAY_EDGES;
}

uint64_t collect_display_edge(int display_bin)
{
        if(display_bin < 0){
                display_bin = 0;
        }
        if(display_bin >= CT_N_DISPLAY_EDGES){
                display_bin = CT_N_DISPLAY_EDGES - 1;
        }
        return display_edges[display_bin];
}

/* Every storage edge from 1000 on is a display edge and below 1000 the display edges are a
 * subset of the storage edges, so the display bin of a storage bin is the last display edge
 * not above the storage bin's lower edge. */
int collect_display_bin(int storage_bin)
{
        uint64_t edge = collect_length_edge(storage_bin);
        int lo = 0;
        int hi = CT_N_DISPLAY_EDGES - 1;

        while(lo < hi){
                int mid = (lo + hi + 1) / 2;
                if(display_edges[mid] <= edge){
                        lo = mid;
                }else{
                        hi = mid - 1;
                }
        }
        return lo;
}

double collect_mean_q(collect_u128 error_sum, uint64_t bases)
{
        double mean_error;

        if(bases == 0){
                return COLLECT_NA;
        }
        mean_error = (double) error_sum / (double) bases / (double) COLLECT_QUAL_ONE;
        if(mean_error <= 0.0){
                return (double) (COLLECT_N_Q - 1);
        }
        return -10.0 * log10(mean_error);
}

/* ------------------------------------------------------------------------------------
 * Collecting (standard-plots.md section 3)
 * ------------------------------------------------------------------------------------ */

#define PHRED_MISSING (-1)

/* Phred value (0..93) of base <pos>, or PHRED_MISSING: QUAL '*' is stored as 0xFF, FASTA
 * has no quality string, and a short quality string covers only part of the read. */
static inline int phred_of(const uint8_t *qs, int qlen, int pos, int offset)
{
        int q;
        if(pos >= qlen){
                return PHRED_MISSING;
        }
        if(qs[pos] == QUAL_MISSING){
                return PHRED_MISSING;
        }
        q = (int) qs[pos] - offset;
        if(q < 0){
                q = 0;
        }
        if(q > COLLECT_N_Q - 1){
                q = COLLECT_N_Q - 1;
        }
        return q;
}

/* Bin of a per-read mean error: the largest q with E[q] >= mean (E never increases with q),
 * integer comparisons only. The mean is at most E[0] = 2^30. */
static inline int mean_q_bin(uint64_t error_sum, uint64_t n)
{
        uint64_t mean = error_sum / n;
        int lo = 0;
        int hi = COLLECT_N_Q - 1;
        while(lo < hi){
                int mid = (lo + hi + 1) / 2;
                if(qual_error[mid] >= mean){
                        lo = mid;
                }else{
                        hi = mid - 1;
                }
        }
        return lo;
}

/* Relative position axis: bin = (100 * i) / L. The reads are walked in increasing position,
 * so the bin is advanced across boundaries (b * L / 100) instead of dividing per base. */
struct rel_walk {
        uint64_t len;   /* aligned length L (> 0) */
        uint64_t thr;   /* (bin + 1) * L: the first 100 * i that belongs to the next bin */
        int bin;
};

static inline void rel_init(struct rel_walk *w, uint64_t len)
{
        w->len = (len > 0) ? len : 1;
        w->thr = w->len;
        w->bin = 0;
}

/* pos must not decrease between calls */
static inline int rel_at(struct rel_walk *w, uint64_t pos)
{
        uint64_t a = 100 * pos;
        while(a >= w->thr){
                w->bin++;
                w->thr += w->len;
        }
        /* a position one past the last base (a deletion at the very end) stays in the last bin */
        return (w->bin < COLLECT_REL_BINS) ? w->bin : COLLECT_REL_BINS - 1;
}

/* The aligned part of a read, both axes, one pass. do_abs and do_aln are constants at the
 * call sites: the first 1000 positions update the absolute axis as well, and reads with an
 * alignment count their bases as aligned. */
struct pos_ctx {
        struct collect_pos_bin *ab;     /* absolute bins of this mate and MAPQ class */
        struct collect_pos_bin *rb;     /* relative bins */
        const uint8_t *seq;
        const uint8_t *qs;
        int qlen;
        int q_offset;
        uint64_t inc_aligned_md;
        uint64_t **legacy_q;            /* legacy Q matrix rows of this MAPQ class */
        int legacy_cols;
        uint64_t rd_error_sum;          /* error sum and count of the whole read (clips too) */
        uint64_t rd_n;
};

static inline __attribute__((always_inline))
void pos_loop(struct pos_ctx *c, struct rel_walk *w, int c5, int from, int to, const int do_abs, const int do_aln)
{
        for(int i = from; i < to; i++){
                int p = c5 + i;
                struct collect_pos_bin *r_ = c->rb + rel_at(w, (uint64_t) i);
                struct collect_pos_bin *a_ = do_abs ? c->ab + i : r_;
                int code = nuc_to_internal[c->seq[p]];
                int q = phred_of(c->qs, c->qlen, p, c->q_offset);

                r_->comp[code]++;
                if(do_aln){
                        r_->aligned++;
                        r_->aligned_md += c->inc_aligned_md;
                }
                if(do_abs){
                        a_->comp[code]++;
                        if(do_aln){
                                a_->aligned++;
                                a_->aligned_md += c->inc_aligned_md;
                        }
                }
                if(q != PHRED_MISSING){
                        uint64_t e = qual_error[q];
                        r_->qual_bases++;
                        r_->qual_error_sum += e;
                        if(do_abs){
                                a_->qual_bases++;
                                a_->qual_error_sum += e;
                                if(i < c->legacy_cols){
                                        c->legacy_q[(q > QUAL_MAX_BIN) ? QUAL_MAX_BIN : q][i]++;
                                }
                        }
                        c->rd_error_sum += e;
                        c->rd_n++;
                }
        }
}

/* Mismatches, insertions and deletions along the alignment strings (read orientation, after
 * fix_orientation). rp is the read position inside the aligned part. An insertion or a
 * deletion is counted once, where the gap opens (2.3 semantics). */
static void collect_alignment(struct stat_collection *s, const struct aln_data *a, int r, int m, int aligned_len)
{
        struct collect_acc *acc = s->acc;
        struct collect_pos_bin *ab = acc->pos_abs[r][m];
        struct collect_pos_bin *rb = acc->pos_rel[r][m];
        struct plot_data *mis_pd = (r == 0) ? s->mis_R1 : s->mis_R2;
        uint64_t **mis = mis_pd->data + COLLECT_N_BASE * m;
        const int mis_cols = mis_pd->len;
        const uint8_t *rd = a->read;
        const uint8_t *gn = a->genome;
        const int aln_len = a->aln_len;
        const int has_mismatch_info = a->has_md_or_x;
        struct rel_walk w;
        int rp = 0;

        rel_init(&w, (uint64_t) aligned_len);
        for(int j = 0; j < aln_len; j++){
                const uint8_t rc = rd[j];
                const uint8_t gc = gn[j];

                if(rc != '-'){
                        if(gc == '-'){
                                if(j == 0 || gn[j - 1] != '-'){
                                        if(rp < COLLECT_ABS_BINS){
                                                ab[rp].insertions++;
                                        }
                                        rb[rel_at(&w, (uint64_t) rp)].insertions++;
                                }
                        }else if(has_mismatch_info && rc != gc){
                                if(rp < COLLECT_ABS_BINS){
                                        ab[rp].mismatches++;
                                }
                                rb[rel_at(&w, (uint64_t) rp)].mismatches++;
                                if(rp < mis_cols){
                                        mis[nuc_to_internal[rc]][rp]++;
                                }
                        }
                        rp++;
                }else if(gc != '-'){
                        if(j == 0 || rd[j - 1] != '-'){
                                if(rp < COLLECT_ABS_BINS){
                                        ab[rp].deletions++;
                                }
                                rb[rel_at(&w, (uint64_t) rp)].deletions++;
                        }
                }
        }
}

/* CIGAR has = or X operations (identity rule 2 of 3.4) */
static int cigar_has_eqx(const struct aln_data *a)
{
        if(a->x_ops > 0){
                return 1;
        }
        for(int i = 0; i < a->n_cigar; i++){
                uint32_t op = bam_cigar_op(a->cigar[i]);
                if(op == BAM_CEQUAL || op == BAM_CDIFF){
                        return 1;
                }
        }
        return 0;
}

/* Gap-compressed identity of a mapped read (3.4) */
static void collect_identity(struct stat_collection *s, const struct aln_data *a, int r, int m)
{
        const uint64_t M = a->m_bases;
        const uint64_t denom = M + a->i_ev + a->d_ev;
        int64_t X = 0;
        int ok = 0;

        if(denom > 0){
                if(a->error >= 0){
                        /* NM counts mismatches and the bases of every gap */
                        X = (int64_t) a->error - (int64_t) a->i_bases - (int64_t) a->d_bases;
                        ok = (X >= 0);
                }else if(cigar_has_eqx(a)){
                        X = (int64_t) a->x_ops;
                        ok = 1;
                }
                if(ok && (uint64_t) X > M){
                        ok = 0;
                }
        }
        if(ok){
                s->acc->identity[r][m][(1000 * (M - (uint64_t) X)) / denom]++;
                s->identity_reads++;
        }else{
                s->identity_na++;
        }
}

/* One primary record. Every decision is made here, per record: quality accumulators iff
 * the QUAL is present, alignment accumulators iff the record has aln_data and is mapped. */
static void collect_record(struct stat_collection *s, struct tl_seq *itm, int q_offset)
{
        struct collect_acc *acc = s->acc;
        struct aln_data *a = (struct aln_data *) itm->data;
        const int len = itm->len;
        const int n = MACRO_MIN(len, (int) itm->seq->len);
        int r = 0;
        int m = MAPQUALBIN_UNMAP;
        int mapped = 0;
        int has_aln = 0;
        int c5 = 0;
        int c3 = 0;
        int bin;

        if(a){
                if(a->flag & BAM_FPAIRED){
                        s->n_paired++;
                        if(a->flag & BAM_FPROPER_PAIR){
                                s->n_proper_paired++;
                        }
                }
                if((a->flag & BAM_FREAD2) && !(a->flag & BAM_FREAD1)){
                        r = 1;
                }
                if(!(a->flag & BAM_FUNMAP)){
                        mapped = 1;
                        m = s->mapq_map->map[a->map_q & 0xFF];
                }
                has_aln = mapped && a->aln_len > 0;
                if(has_aln){
                        c5 = (int) a->n_clip5;
                        c3 = (int) a->n_clip3;
                        if(c5 < 0 || c3 < 0 || c5 + c3 > n){
                                c5 = 0;
                                c3 = 0;
                        }
                }
        }

        s->basic_nums[r][m]++;

        /* length (3.1): the full read, clipped bases included */
        s->total_bases += (uint64_t) len;
        if((uint64_t) len < s->min_length){
                s->min_length = (uint64_t) len;
        }
        if((uint64_t) len > s->max_length){
                s->max_length = (uint64_t) len;
        }
        bin = collect_length_bin((uint64_t) len);
        acc->len_reads[r][m][bin]++;
        acc->len_bases[r][m][bin] += (uint64_t) len;
        {
                /* legacy length histogram: clip-trimmed, longer reads in the last column */
                struct plot_data *ld = (r == 0) ? s->len_dist_R1 : s->len_dist_R2;
                ld->data[m][clamp_idx(len - (c5 + c3), ld->len)]++;
        }

        if(mapped){
                s->mapped_bases += (uint64_t) len;
                s->clipped_bases += a->s_bases;
        }

        /* bases: clipped bases only feed the read mean Q, the aligned part feeds everything */
        {
                const int aligned_len = n - c5 - c3;
                struct pos_ctx c;
                struct rel_walk w;

                c.ab = acc->pos_abs[r][m];
                c.rb = acc->pos_rel[r][m];
                c.seq = (const uint8_t *) itm->seq->str;
                c.qs = (const uint8_t *) itm->qual->str;
                c.qlen = (int) itm->qual->len;
                c.q_offset = q_offset;
                c.inc_aligned_md = (has_aln && a->has_md_or_x) ? 1 : 0;
                {
                        struct plot_data *lq = (r == 0) ? s->qual_comp_R1 : s->qual_comp_R2;
                        c.legacy_q = lq->data + (lq->L / COLLECT_N_MAPQ) * m;
                        c.legacy_cols = lq->len;
                }
                c.rd_error_sum = 0;
                c.rd_n = 0;

                for(int p = 0; p < c5; p++){
                        int q = phred_of(c.qs, c.qlen, p, q_offset);
                        if(q != PHRED_MISSING){
                                c.rd_error_sum += qual_error[q];
                                c.rd_n++;
                        }
                }
                for(int p = n - c3; p < n; p++){
                        int q = phred_of(c.qs, c.qlen, p, q_offset);
                        if(q != PHRED_MISSING){
                                c.rd_error_sum += qual_error[q];
                                c.rd_n++;
                        }
                }
                if(aligned_len > 0){
                        const int n_abs = MACRO_MIN(aligned_len, COLLECT_ABS_BINS);
                        rel_init(&w, (uint64_t) aligned_len);
                        if(has_aln){
                                pos_loop(&c, &w, c5, 0, n_abs, 1, 1);
                                pos_loop(&c, &w, c5, n_abs, aligned_len, 0, 1);
                        }else{
                                pos_loop(&c, &w, c5, 0, n_abs, 1, 0);
                                pos_loop(&c, &w, c5, n_abs, aligned_len, 0, 0);
                        }
                        if(aligned_len > COLLECT_ABS_BINS){
                                s->bases_beyond_position_axis += (uint64_t) (aligned_len - COLLECT_ABS_BINS);
                        }
                }
                if(c.rd_n){
                        acc->read_q[r][m][mean_q_bin(c.rd_error_sum, c.rd_n)]++;
                        s->qual_bases_all += c.rd_n;
                        s->qual_error_sum_all += c.rd_error_sum;
                }
                if(has_aln){
                        collect_alignment(s, a, r, m, aligned_len);
                        if(!a->has_md_or_x){
                                s->reads_without_md++;
                        }
                }
        }
        if(mapped){
                collect_identity(s, a, r, m);
        }
}

int collect_stats(struct tl_seq_buffer *sb, struct stat_collection *s)
{
        int max_len = 0;
        int pos_cols;
        int len_cols;

        for(int i = 0; i < sb->num_seq; i++){
                if(sb->sequences[i]->len > max_len){
                        max_len = sb->sequences[i]->len;
                }
        }
        /* legacy matrices: as wide as the longest read, but never wider than the cap */
        pos_cols = MACRO_MIN(max_len + 1, COLLECT_LEGACY_POS_COLS);
        len_cols = MACRO_MIN(max_len + 1, COLLECT_LEGACY_LEN_COLS);

        RUN(plot_data_resize_len(s->base_comp_R1, pos_cols));
        RUN(plot_data_resize_len(s->base_comp_R2, pos_cols));

        RUN(plot_data_resize_len(s->qual_comp_R1, pos_cols));
        RUN(plot_data_resize_len(s->qual_comp_R2, pos_cols));

        RUN(plot_data_resize_len(s->mis_R1, pos_cols));
        RUN(plot_data_resize_len(s->mis_R2, pos_cols));

        RUN(plot_data_resize_len(s->ins_R1, pos_cols));
        RUN(plot_data_resize_len(s->ins_R2, pos_cols));

        RUN(plot_data_resize_len(s->del_R1, pos_cols));
        RUN(plot_data_resize_len(s->del_R2, pos_cols));

        RUN(plot_data_resize_len(s->len_dist_R1, len_cols));
        RUN(plot_data_resize_len(s->len_dist_R2, len_cols));

        if(s->collect_end){
                RUN(plot_data_resize_len(s->base_comp_R1_end, pos_cols));
                RUN(plot_data_resize_len(s->base_comp_R2_end, pos_cols));

                RUN(plot_data_resize_len(s->qual_comp_R1_end, pos_cols));
                RUN(plot_data_resize_len(s->qual_comp_R2_end, pos_cols));
        }

        for(int i = 0; i < sb->num_seq; i++){
                collect_record(s, sb->sequences[i], sb->base_quality_offset);
        }

        /* optional legacy plots (unreachable since --plotend and --norm set nothing) */
        if(s->collect_end){
                RUN(collect_stats_end_seqqual(sb, s));
        }
        if(s->norm_len_plot){
                RUN(collect_stats_norm_seqqual(sb, s));
        }
        return OK;
ERROR:
        return FAIL;
}

int collect_stats_norm_seqqual(struct tl_seq_buffer *sb, struct stat_collection *s)
{
        uint64_t** t= NULL;
        uint64_t** t_q= NULL;
        int32_t q_idx;

        q_idx = MAPQUALBIN_UNMAP;
        for(int i = 0; i < sb->num_seq;i++){
                uint8_t* seq = NULL;
                struct tl_seq *itm = sb->sequences[i];
                int clip_start = 0;
                int clip_end = 0;
                q_idx = MAPQUALBIN_UNMAP;
                int read = 1;
                if(itm->data){
                        struct aln_data* a = NULL;

                        a = itm->data;
                        clip_start = a->n_clip5;
                        clip_end = a->n_clip3;

                        q_idx = s->mapq_map->map[ a->map_q];
                        if(a->flag  & BAM_FREAD1){
                                read = 1;
                        }else if(a->flag  & BAM_FREAD2){
                                read = 2;
                        }
                        if(a->flag & BAM_FUNMAP){
                                q_idx = MAPQUALBIN_UNMAP;
                        }
                }
                seq = itm->seq->str;

                /* start collecting...  */
                if(read == 1){
                        t = s->base_comp_len_norm_R1->data  + s->base_comp_len_norm_R1->group_size * q_idx;
                        t_q = s->qual_comp_len_norm_R1->data  + s->qual_comp_len_norm_R1->group_size * q_idx;
                }else{
                        t = s->base_comp_len_norm_R2->data  + s->base_comp_len_norm_R2->group_size * q_idx;
                        t_q = s->qual_comp_len_norm_R2->data  + s->qual_comp_len_norm_R2->group_size * q_idx;
                }

                /* ASSERT(s->base_comp_R1_end->len ==s->qual_comp_R1_end->len,"Lengths differ seq / qual"); */
                double plot_len = (double) s->base_comp_len_norm_R1->len;


                int len = MACRO_MIN(itm->len - clip_end, itm->seq->len);
                int t_len = itm->len - (clip_start + clip_end);
                int p = 0;
                if(t_len <= 0){
                        continue;
                }
                for(int j = MACRO_MAX(clip_start, 0);j < len;j++){
                        int idx = (int) floor((double) p / (double) (t_len) * plot_len);
                        /* rounding can push the last positions onto plot_len */
                        idx = clamp_idx(idx, (int) plot_len);
                        /* LOG_MSG("%d %d -> idx: %d", j,p,idx); */
                        t[nuc_to_internal[seq[j]]][idx]++;

                        int q = get_base_qual(itm, j, sb->base_quality_offset);
                        if(q >= 0){
                                t_q[q][idx]++;
                        }
                        p++;
                }
        }
        return OK;
}


int collect_stats_end_seqqual(struct tl_seq_buffer *sb, struct stat_collection *s)
{
        uint64_t** t= NULL;
        uint64_t** t_q= NULL;
        int32_t q_idx;

        q_idx = MAPQUALBIN_UNMAP;
        for(int i = 0; i < sb->num_seq;i++){
                uint8_t* seq = NULL;
                struct tl_seq *itm = sb->sequences[i];

                q_idx = MAPQUALBIN_UNMAP;
                int read = 1;
                if(itm->data){
                        struct aln_data* a = NULL;
                        a = itm->data;
                        q_idx = s->mapq_map->map[ a->map_q];
                        if(a->flag  & BAM_FREAD1){
                                read = 1;
                        }else if(a->flag  & BAM_FREAD2){
                                read = 2;
                        }
                        if(a->flag & BAM_FUNMAP){
                                q_idx = MAPQUALBIN_UNMAP;
                        }
                }
                seq = itm->seq->str;

                /* start collecting...  */
                if(read == 1){
                        t = s->base_comp_R1_end->data  + s->base_comp_R1->group_size * q_idx;
                        t_q = s->qual_comp_R1_end->data  + s->qual_comp_R1->group_size * q_idx;
                }else{
                        t = s->base_comp_R2_end->data  + s->base_comp_R2->group_size * q_idx;
                        t_q = s->qual_comp_R2_end->data  + s->qual_comp_R2->group_size * q_idx;
                }

                ASSERT(s->base_comp_R1_end->len ==s->qual_comp_R1_end->len,"Lengths differ seq / qual");
                int plot_len = s->base_comp_R1_end->len;
                int len = MACRO_MIN(itm->len, itm->seq->len) - 1;
                int idx = 0;

                while(len >= 0 && idx < plot_len){
                        /* LOG_MSG("idx = %d", idx); */
                        t[nuc_to_internal[seq[len]]][idx]++;
                        int q = get_base_qual(itm, len, sb->base_quality_offset);
                        if(q >= 0){
                                t_q[q][idx]++;
                        }
                        idx++;
                        len--;
                }
        }
        return OK;
ERROR:
        return FAIL;
}



/* ------------------------------------------------------------------------------------
 * Merge (standard-plots.md 3.6)
 * ------------------------------------------------------------------------------------ */

/* every legacy matrix of a stat_collection */
#define LEGACY_PLOTS(X)                                                 \
        X(base_comp_R1) X(base_comp_R2) X(qual_comp_R1) X(qual_comp_R2) \
        X(mis_R1) X(mis_R2) X(ins_R1) X(ins_R2) X(del_R1) X(del_R2)     \
        X(len_dist_R1) X(len_dist_R2)                                   \
        X(base_comp_R1_end) X(base_comp_R2_end)                         \
        X(qual_comp_R1_end) X(qual_comp_R2_end)                         \
        X(base_comp_len_norm_R1) X(qual_comp_len_norm_R1)               \
        X(base_comp_len_norm_R2) X(qual_comp_len_norm_R2)

/* Add the counts of <src> to <dst>. The matrices of the two sides may differ in width (each
 * is as wide as the longest read it saw): <dst> grows to the wider one. */
static int merge_plot_data(struct plot_data *dst, const struct plot_data *src)
{
        if(src == NULL){
                return OK;
        }
        ASSERT(dst != NULL, "Cannot merge a legacy matrix into a missing one");
        ASSERT(dst->L == src->L, "Legacy matrices differ in height");
        RUN(plot_data_resize_len(dst, src->len));
        for(int i = 0; i < src->L; i++){
                for(int j = 0; j < src->len; j++){
                        dst->data[i][j] += src->data[i][j];
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int stat_collection_merge(struct stat_collection *dst, const struct stat_collection *src)
{
        uint64_t *d = NULL;
        const uint64_t *x = NULL;

        ASSERT(dst != NULL && src != NULL, "No stat collection");
        ASSERT(dst != src, "Cannot merge a collection into itself");

        /* struct collect_acc holds nothing but uint64_t counters */
        _Static_assert(sizeof(struct collect_acc) % sizeof(uint64_t) == 0, "collect_acc is all uint64_t");
        d = (uint64_t *) dst->acc;
        x = (const uint64_t *) src->acc;
        for(size_t i = 0; i < sizeof(struct collect_acc) / sizeof(uint64_t); i++){
                d[i] += x[i];
        }

        for(int r = 0; r < COLLECT_N_MATE; r++){
                for(int m = 0; m < COLLECT_N_MAPQ; m++){
                        dst->basic_nums[r][m] += src->basic_nums[r][m];
                }
        }

        dst->n_total += src->n_total;
        dst->n_read1 += src->n_read1;
        dst->n_read2 += src->n_read2;
        dst->n_paired += src->n_paired;
        dst->n_proper_paired += src->n_proper_paired;
        dst->n_secondary += src->n_secondary;
        dst->n_supplementary += src->n_supplementary;

        dst->total_bases += src->total_bases;
        if(src->min_length < dst->min_length){
                dst->min_length = src->min_length;
        }
        if(src->max_length > dst->max_length){
                dst->max_length = src->max_length;
        }
        dst->mapped_bases += src->mapped_bases;
        dst->clipped_bases += src->clipped_bases;
        dst->qual_bases_all += src->qual_bases_all;
        dst->qual_error_sum_all += src->qual_error_sum_all;
        dst->reads_without_md += src->reads_without_md;
        dst->identity_reads += src->identity_reads;
        dst->identity_na += src->identity_na;
        dst->bases_beyond_position_axis += src->bases_beyond_position_axis;
        dst->has_sq |= src->has_sq;
        dst->is_partial_report |= src->is_partial_report;

#define MERGE_PLOT(NAME) RUN(merge_plot_data(dst->NAME, src->NAME));
        LEGACY_PLOTS(MERGE_PLOT)
#undef MERGE_PLOT
        return OK;
ERROR:
        return FAIL;
}

/* ------------------------------------------------------------------------------------
 * Finalise (standard-plots.md 3.7)
 * ------------------------------------------------------------------------------------ */

/* The legacy base composition, insertion and deletion matrices hold exactly the first
 * 1000 positions of the absolute axis, so they are copied from there instead of being
 * counted a second time per base. */
static void derive_legacy_matrices(struct stat_collection *s)
{
        for(int r = 0; r < COLLECT_N_MATE; r++){
                struct plot_data *bc = (r == 0) ? s->base_comp_R1 : s->base_comp_R2;
                struct plot_data *ins = (r == 0) ? s->ins_R1 : s->ins_R2;
                struct plot_data *del = (r == 0) ? s->del_R1 : s->del_R2;
                const int cols = MACRO_MIN(MACRO_MIN(bc->len, MACRO_MIN(ins->len, del->len)), COLLECT_ABS_BINS);

                for(int m = 0; m < COLLECT_N_MAPQ; m++){
                        const struct collect_pos_bin *pb = s->acc->pos_abs[r][m];
                        for(int i = 0; i < cols; i++){
                                for(int b = 0; b < COLLECT_N_BASE; b++){
                                        bc->data[COLLECT_N_BASE * m + b][i] = pb[i].comp[b];
                                }
                                ins->data[m][i] = pb[i].insertions;
                                del->data[m][i] = pb[i].deletions;
                        }
                }
        }
}

int stat_collection_finalise(struct stat_collection *s)
{
        uint64_t n_unmapped;
        uint64_t n_mapped;

        ASSERT(s != NULL,"No stat collection");

        s->n_total = 0;
        s->n_read1 = 0;
        s->n_read2 = 0;
        for(int i = 0; i < COLLECT_N_MATE; i++){
                for(int j = 0;j < s->mapq_map->n_bin;j++){
                        s->n_total += s->basic_nums[i][j];
                        if(i == 0){
                                s->n_read1 += s->basic_nums[i][j];
                        }else{
                                s->n_read2 += s->basic_nums[i][j];
                        }
                }
        }

        if(s->n_read2 && !s->titles_done){
                RUN(tld_append(s->base_comp_R1->title," Read 1"));;
                RUN(tld_append(s->base_comp_R2->title," Read 2"));

                RUN(tld_append(s->qual_comp_R1->title," Read 1"));;
                RUN(tld_append(s->qual_comp_R2->title," Read 2"));

                RUN(tld_append(s->mis_R1->title," Read 1"));;
                RUN(tld_append(s->mis_R2->title," Read 2"));

                RUN(tld_append(s->ins_R1->title," Read 1"));;
                RUN(tld_append(s->ins_R2->title," Read 2"));

                RUN(tld_append(s->del_R1->title," Read 1"));;
                RUN(tld_append(s->del_R2->title," Read 2"));

                RUN(tld_append(s->len_dist_R1->title," Read 1"));;
                RUN(tld_append(s->len_dist_R2->title," Read 2"));
        }
        s->titles_done = 1;

        derive_legacy_matrices(s);

        if(s->n_total == 0){
                /* no read: min_length still holds its UINT64_MAX start value */
                s->min_length = 0;
        }

        n_unmapped = s->basic_nums[0][MAPQUALBIN_UNMAP] + s->basic_nums[1][MAPQUALBIN_UNMAP];
        n_mapped = s->n_total - n_unmapped;

        /* is_aligned = has_sq || at least one primary record is mapped (3.6) */
        s->is_aligned = (s->has_sq || n_mapped > 0) ? 1 : 0;

        /* length: N50 from the longest end, median from the shortest end, both as the lower
         * edge of the storage bin (exact below 1000 bp) */
        {
                uint64_t reads_at[COLLECT_N_LEN_BINS];
                uint64_t bases_at[COLLECT_N_LEN_BINS];
                uint64_t cum = 0;

                for(int b = 0; b < COLLECT_N_LEN_BINS; b++){
                        reads_at[b] = 0;
                        bases_at[b] = 0;
                        for(int r = 0; r < COLLECT_N_MATE; r++){
                                for(int m = 0; m < COLLECT_N_MAPQ; m++){
                                        reads_at[b] += s->acc->len_reads[r][m][b];
                                        bases_at[b] += s->acc->len_bases[r][m][b];
                                }
                        }
                }
                s->median_length = 0;
                s->n50 = 0;
                if(s->n_total){
                        for(int b = 0; b < COLLECT_N_LEN_BINS; b++){
                                cum += reads_at[b];
                                if(2 * cum >= s->n_total){
                                        s->median_length = collect_length_edge(b);
                                        break;
                                }
                        }
                }
                if(s->total_bases){
                        cum = 0;
                        for(int b = COLLECT_N_LEN_BINS - 1; b >= 0; b--){
                                cum += bases_at[b];
                                if(2 * cum >= s->total_bases){
                                        s->n50 = collect_length_edge(b);
                                        break;
                                }
                        }
                }
        }
        s->layout = (s->median_length >= COLLECT_LAYOUT_MIN_LONG_MEDIAN) ? COLLECT_LAYOUT_LONG : COLLECT_LAYOUT_SHORT;

        s->mean_q = collect_mean_q(s->qual_error_sum_all, s->qual_bases_all);

        s->pct_mapped = COLLECT_NA;
        s->pct_mapq_ge30 = COLLECT_NA;
        if(s->is_aligned && s->n_total){
                uint64_t n_ge30 = s->basic_nums[0][MAPQUALBIN_gt30] + s->basic_nums[1][MAPQUALBIN_gt30];
                s->pct_mapped = 100.0 * (double) n_mapped / (double) s->n_total;
                s->pct_mapq_ge30 = 100.0 * (double) n_ge30 / (double) s->n_total;
        }
        s->pct_clipped = COLLECT_NA;
        if(s->mapped_bases){
                s->pct_clipped = 100.0 * (double) s->clipped_bases / (double) s->mapped_bases;
        }

        /* median identity: lower edge of the bin holding the lower median read */
        s->median_identity = COLLECT_NA;
        if(s->identity_reads){
                uint64_t cum = 0;
                for(int b = 0; b < COLLECT_N_IDENT; b++){
                        for(int r = 0; r < COLLECT_N_MATE; r++){
                                for(int m = 0; m < COLLECT_N_MAPQ; m++){
                                        cum += s->acc->identity[r][m][b];
                                }
                        }
                        if(2 * cum >= s->identity_reads){
                                s->median_identity = (double) b / 10.0;
                                break;
                        }
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int stat_collection_config_additional_plot(struct stat_collection *s,struct samstat_param *p)
{
        char** base_label = NULL;

        galloc(&base_label, 5,2);
        snprintf(base_label[0], 2, "%s","A");
        snprintf(base_label[1], 2, "%s","C");
        snprintf(base_label[2], 2, "%s","G");
        snprintf(base_label[3], 2, "%s","T");
        snprintf(base_label[4], 2, "%s","N");

        if(p->norm_len_plot){
                s->norm_len_plot = 1;
                /* how should I determine the number + size of bins? */

                /* I need base and base distripbution plot  */
                RUN(plot_data_alloc(&s->base_comp_len_norm_R1, 100, 5* s->mapq_map->n_bin));
                RUN(plot_data_alloc(&s->qual_comp_len_norm_R1, 100, 42* s->mapq_map->n_bin));

                RUN(plot_data_alloc(&s->base_comp_len_norm_R2, 100, 5* s->mapq_map->n_bin));
                RUN(plot_data_alloc(&s->qual_comp_len_norm_R2, 100, 42* s->mapq_map->n_bin));


                /* Need to add config */
                RUN(plot_data_config(s->base_comp_len_norm_R1,
                                     PLOT_TYPE_BAR,
                                     PLOT_MOD_NORMAL,
                                     5,
                                     PLOT_EMPTY_SERIES, /* plot empty series - yes please  */
                                     "bc_R1_norm",
                                     "Base composition",
                                     "Base composition split into 1% bins along the read length.",
                                     "Counts",
                                     "BaseCompositionnorm",
                                     base_label,
                                     s->mapq_map->description
                            ));

                RUN(plot_data_config(s->base_comp_len_norm_R2,
                                     PLOT_TYPE_BAR,
                                     PLOT_MOD_NORMAL,
                                     5,
                                     PLOT_EMPTY_SERIES, /* plot empty series - yes please  */
                                     "bc_R2_norm",
                                     "Base composition",
                                     "Base composition split into 1% bins along the read length.",
                                     "Counts",
                                     "BaseCompositionnorm",
                                     base_label,
                                     s->mapq_map->description
                            ));

                RUN(plot_data_config(s->qual_comp_len_norm_R1,
                                     PLOT_TYPE_LINES,
                                     PLOT_MOD_ERROR_BAR,
                                     42,
                                     0,
                                     "bq_norm_R1",
                                     "Base quality distribution",
                                     "Base quality distribution split into 1% bins along the read length.",
                                     "Base Quality",
                                     "QualCompositionnorm",
                                     NULL,
                                     s->mapq_map->description
                            ));

                RUN(plot_data_config(s->qual_comp_len_norm_R2,
                                     PLOT_TYPE_LINES,
                                     PLOT_MOD_ERROR_BAR,
                                     42,
                                     0,
                                     "bq_norm_R2",
                                     "Base quality distribution",
                                     "Base quality distribution split into 1% bins along the read length.",
                                     "Base Quality",
                                     "QualCompositionR2norm",
                                     NULL,
                                     s->mapq_map->description
                            ));


        }
        if(p->collect_end){
                s->collect_end = 1;
                RUN(plot_data_alloc(&s->base_comp_R1_end, 250,  5* s->mapq_map->n_bin));
                RUN(plot_data_alloc(&s->base_comp_R2_end, 250,  5* s->mapq_map->n_bin));

                /* Alloc 250 len and 42 slots (base qualities ) for every mapq bin  */
                RUN(plot_data_alloc(&s->qual_comp_R1_end, 250,  42* s->mapq_map->n_bin));
                RUN(plot_data_alloc(&s->qual_comp_R2_end, 250,  42* s->mapq_map->n_bin));



                /* s->mapq_map->description */
                RUN(plot_data_config(s->base_comp_R1_end,
                                     PLOT_TYPE_BAR,
                                     PLOT_MOD_NORMAL,
                                     5,
                                     PLOT_EMPTY_SERIES, /* plot empty series - yes please  */
                                     "bc_R1_end",
                                     "Base Composition at read end",
                                     "Position relative to read end",
                                     "Counts",
                                     "BaseCompositionend",
                                     base_label,
                                     s->mapq_map->description
                            ));

                RUN(plot_data_config(s->base_comp_R2_end,
                                     PLOT_TYPE_BAR,
                                     PLOT_MOD_NORMAL,
                                     5,
                                     PLOT_EMPTY_SERIES, /* plot empty series - yes please  */
                                     "bc_R2_end",
                                     "Base Composition at read end",
                                     "Position relative to read end",
                                     "Counts",
                                     "BaseCompositionR2end",
                                     base_label,
                                     s->mapq_map->description
                            ));



                /* qual composition */

                RUN(plot_data_config(s->qual_comp_R1_end,
                                     PLOT_TYPE_LINES,
                                     PLOT_MOD_ERROR_BAR,
                                     42,
                                     0,
                                     "bq_R1_end",
                                     "Base quality distribution at read end",
                                     "Position relative to read end",
                                     "Base Quality",
                                     "QualCompositionend",
                                     NULL,
                                     s->mapq_map->description
                            ));

                RUN(plot_data_config(s->qual_comp_R2_end ,
                                     PLOT_TYPE_LINES,
                                     PLOT_MOD_ERROR_BAR,
                                     42,
                                     0,
                                     "bq_R2_end",
                                     "Base quality distribution at read end",
                                     "Position relative to read end",
                                     "Base Quality",
                                     "QualCompositionR2end",
                                     NULL,
                                     s->mapq_map->description
                            ));

        }
        gfree(base_label);
        return OK;
ERROR:
        gfree(base_label);
        return FAIL;
}


int stat_collection_alloc(struct stat_collection **stats)
{
        struct stat_collection* s = NULL;

        MMALLOC(s,sizeof(struct stat_collection)) ;

        s->base_comp_R1 = NULL;
        s->base_comp_R2 = NULL;

        s->qual_comp_R1 = NULL;
        s->qual_comp_R2 = NULL;

        s->mis_R1 = NULL;
        s->mis_R2 = NULL;

        s->ins_R1 = NULL;
        s->ins_R2 = NULL;

        s->del_R1 = NULL;
        s->del_R2 = NULL;

        s->len_dist_R1 = NULL;
        s->len_dist_R2 = NULL;

        /* optional plots  */
        s->base_comp_R1_end = NULL;
        s->base_comp_R2_end = NULL;

        s->qual_comp_R1_end = NULL;
        s->qual_comp_R2_end = NULL;
        /* norm_len_plot */

        s->base_comp_len_norm_R1 = NULL;
        s->qual_comp_len_norm_R1 = NULL;

        s->base_comp_len_norm_R2 = NULL;
        s->qual_comp_len_norm_R2 = NULL;


        s->mapq_map = NULL;
        s->acc = NULL;

        s->basic_nums = NULL;

        s->n_paired = 0;
        s->n_proper_paired = 0;
        s->n_secondary = 0;
        s->n_supplementary = 0;

        s->n_read1 = 0;
        s->n_read2 = 0;

        s->is_aligned = 0;
        s->is_partial_report = 0;
        s->collect_end = 0;

        /* v2 accumulators: zero filled by MMALLOC */
        s->acc = NULL;
        MMALLOC(s->acc, sizeof(struct collect_acc));
        s->min_length = UINT64_MAX;
        s->mean_q = COLLECT_NA;
        s->pct_mapped = COLLECT_NA;
        s->pct_mapq_ge30 = COLLECT_NA;
        s->pct_clipped = COLLECT_NA;
        s->median_identity = COLLECT_NA;
        RUN(get_mapqual_bins(&s->mapq_map));


        galloc(&s->basic_nums, 2, s->mapq_map->n_bin);
        for(int i = 0; i < 2;i++){
                for(int j = 0; j < s->mapq_map->n_bin;j++){
                        s->basic_nums[i][j] = 0;
                }
        }

        /* Alloc 250 len and 5 slots (ACTGN) for every mapq bin  */
        RUN(plot_data_alloc(&s->base_comp_R1, 250,  5* s->mapq_map->n_bin));
        RUN(plot_data_alloc(&s->base_comp_R2, 250,  5* s->mapq_map->n_bin));

        /* Alloc 250 len and 42 slots (base qualities ) for every mapq bin  */
        RUN(plot_data_alloc(&s->qual_comp_R1, 250,  42* s->mapq_map->n_bin));
        RUN(plot_data_alloc(&s->qual_comp_R2, 250,  42* s->mapq_map->n_bin));

        /* Alloc 250 len and 5 slots (ACTGN) for every mapq bin  */
        RUN(plot_data_alloc(&s->mis_R1, 250,  5* s->mapq_map->n_bin));
        RUN(plot_data_alloc(&s->mis_R2, 250,  5* s->mapq_map->n_bin));

        /* Alloc 250 len and 1 slot for every mapq bin */
        RUN(plot_data_alloc(&s->ins_R1, 250,  s->mapq_map->n_bin));
        RUN(plot_data_alloc(&s->ins_R2, 250,  s->mapq_map->n_bin));

        /* Alloc 250 len and 1 slot for every mapq bin */
        RUN(plot_data_alloc(&s->del_R1, 250,  s->mapq_map->n_bin));
        RUN(plot_data_alloc(&s->del_R2, 250,  s->mapq_map->n_bin));

        /* Alloc 250 len and 1 slot for every mapq bin */
        RUN(plot_data_alloc(&s->len_dist_R1, 250,  s->mapq_map->n_bin));
        RUN(plot_data_alloc(&s->len_dist_R2, 250,  s->mapq_map->n_bin));

        /* Configure plots  */

        /* base composition */

        char** base_label = NULL;

        galloc(&base_label, 5,2);
        snprintf(base_label[0], 2, "%s","A");
        snprintf(base_label[1], 2, "%s","C");
        snprintf(base_label[2], 2, "%s","G");
        snprintf(base_label[3], 2, "%s","T");
        snprintf(base_label[4], 2, "%s","N");


        /* s->mapq_map->description */
        RUN(plot_data_config(s->base_comp_R1,
                             PLOT_TYPE_BAR,
                             PLOT_MOD_NORMAL,
                             5,
                             PLOT_EMPTY_SERIES, /* plot empty series - yes please  */
                             "bc_R1",
                             "Base Composition",
                             "Length",
                             "Counts",
                             "BaseComposition",
                             base_label,
                             s->mapq_map->description
                    ));

        RUN(plot_data_config(s->base_comp_R2,
                             PLOT_TYPE_BAR,
                             PLOT_MOD_NORMAL,
                             5,
                             PLOT_EMPTY_SERIES, /* plot empty series - yes please  */
                             "bc_R2",
                             "Base Composition",
                             "Length",
                             "Counts",
                             "BaseCompositionR2",
                             base_label,
                             s->mapq_map->description
                    ));



        /* qual composition */

        RUN(plot_data_config(s->qual_comp_R1,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_ERROR_BAR,
                             42,
                             0,
                             "bq_R1",
                             "Base quality distribution",
                             "Length",
                             "Base Quality",
                             "QualComposition",
                             NULL,
                             s->mapq_map->description
                    ));

        RUN(plot_data_config(s->qual_comp_R2 ,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_ERROR_BAR,
                             42,
                             0,
                             "bq_R2",
                             "Base quality distribution",
                             "Length",
                             "Base Quality",
                             "QualCompositionR2",
                             NULL,
                             s->mapq_map->description
                    ));


        /* mismatch composition */
        /* LOG_MSG("MIS comp"); */
        RUN(plot_data_config(s->mis_R1,
                             PLOT_TYPE_BAR,
                             PLOT_MOD_NORMAL,
                             5,
                             PLOT_EMPTY_SERIES,
                             "mis_R1",
                             "Distribution of Mismatches",
                             "Length",
                             "Counts",
                             "MismatchDistribution",
                             base_label,
                             s->mapq_map->description
                    ));


        RUN(plot_data_config(s->mis_R2,
                             PLOT_TYPE_BAR,
                             PLOT_MOD_NORMAL,
                             5,
                             PLOT_EMPTY_SERIES,
                             "mis_R2",
                             "Distribution of Mismatches",
                             "Length",
                             "Counts",
                             "MismatchDistributionR2",
                             base_label,
                             s->mapq_map->description
                    ));




        /* ins composition */
        /* LOG_MSG("Ins comp"); */
        RUN(plot_data_config(s->ins_R1,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_NORMAL,
                             0,//s->mapq_map->n_bin,
                             0,
                             "ins_R1",
                             "Distribution of Insertions",
                             "Length",
                             "Counts",
                             "InsertionDistribution",
                             s->mapq_map->description,
                             NULL
                    ));

        RUN(plot_data_config(s->ins_R2,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_NORMAL,
                             0,//s->mapq_map->n_bin,
                             0,
                             "ins_R2",
                             "Distribution of Insertions",
                             "Length",
                             "Counts",
                             "InsertionDistribution",
                             s->mapq_map->description,
                             NULL
                    ));

        /* del composition */
        /* LOG_MSG("Del comp"); */
        RUN(plot_data_config(s->del_R1,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_NORMAL,
                             0,//s->mapq_map->n_bin,
                             0,
                             "del_R1",
                             "Distribution of Deletions",
                             "Length",
                             "Counts",
                             "DeletionDistribution",
                             s->mapq_map->description,
                             NULL
                    ));

        RUN(plot_data_config(s->del_R2,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_NORMAL,
                             0,//s->mapq_map->n_bin,
                             0,
                             "del_R2",
                             "Distribution of Deletions",
                             "Length",
                             "Counts",
                             "DeletionDistribution",
                             s->mapq_map->description,
                             NULL
                    ));

        /* LOG_MSG("LEN IDST "); */
        /* Length dist  */
        RUN(plot_data_config(s->len_dist_R1,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_DENSITY,
                             0,//s->mapq_map->n_bin,
                             0,
                             "len_R1",
                             "Length Distribution",
                             "Length",
                             "Counts",
                             "LengthDistribution",
                             s->mapq_map->description,
                             NULL

                    ));

        RUN(plot_data_config(s->len_dist_R2,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_DENSITY,
                             0,//s->mapq_map->n_bin,
                             0,
                             "len_R2",
                             "Length Distribution",
                             "Length",
                             "Counts",
                             "LengthDistribution",
                             s->mapq_map->description,
                             NULL
                    ));

        /* fprintf(stderr,"%s -> %d\n", TLD_STR(s->len_dist_R1->title), s->len_dist_R1->group_size); */
        /* fprintf(stderr,"%s -> %d\n", TLD_STR(s->base_comp_R1->title), s->base_comp_R1->group_size); */
        /* fprintf(stderr,"%s -> %d\n", TLD_STR(s->qual_comp_R1->title), s->qual_comp_R1->group_size); */
        /* fprintf(stderr,"%s -> %d\n", TLD_STR(s->mis_R1->title), s->mis_R1->group_size); */
        /* fprintf(stderr,"%s -> %d\n", TLD_STR(s->ins_R1->title), s->ins_R1->group_size); */
        /* fprintf(stderr,"%s -> %d\n", TLD_STR(s->del_R1->title), s->del_R1->group_size); */


        /* exit(0); */
        gfree(base_label);
        *stats = s;

        return OK;
ERROR:
        return FAIL;
}

void stat_collection_free(struct stat_collection *s)
{
        if(s){

                if(s->qual_comp_len_norm_R1){
                        plot_data_free(s->qual_comp_len_norm_R1);
                }
                if(s->base_comp_len_norm_R1){
                        plot_data_free(s->base_comp_len_norm_R1);
                }
                if(s->qual_comp_len_norm_R2){
                        plot_data_free(s->qual_comp_len_norm_R2);
                }
                if(s->base_comp_len_norm_R2){
                        plot_data_free(s->base_comp_len_norm_R2);
                }


                if(s->basic_nums){
                        gfree(s->basic_nums);
                }
                if(s->base_comp_R1){
                        plot_data_free(s->base_comp_R1);
                }
                if(s->base_comp_R2){
                        plot_data_free(s->base_comp_R2);
                }

                if(s->qual_comp_R1){
                        plot_data_free(s->qual_comp_R1);
                }
                if(s->qual_comp_R2){
                        plot_data_free(s->qual_comp_R2);
                }

                if(s->mis_R1){
                        plot_data_free(s->mis_R1);
                }
                if(s->mis_R2){
                        plot_data_free(s->mis_R2);
                }

                if(s->ins_R1){
                        plot_data_free(s->ins_R1);
                }
                if(s->ins_R2){
                        plot_data_free(s->ins_R2);
                }

                if(s->del_R1){
                        plot_data_free(s->del_R1);
                }
                if(s->del_R2){
                        plot_data_free(s->del_R2);
                }

                if(s->len_dist_R1){
                        plot_data_free(s->len_dist_R1);
                }
                if(s->len_dist_R2){
                        plot_data_free(s->len_dist_R2);
                }
                /* Additional plots */

                if(s->base_comp_R1_end){
                        plot_data_free(s->base_comp_R1_end);
                }
                if(s->base_comp_R2_end){
                        plot_data_free(s->base_comp_R2_end);
                }

                if(s->qual_comp_R1_end){
                        plot_data_free(s->qual_comp_R1_end);
                }
                if(s->qual_comp_R2_end){
                        plot_data_free(s->qual_comp_R2_end);
                }

                if(s->mapq_map){
                        free_mapqual_bins(s->mapq_map);
                }
                if(s->acc){
                        MFREE(s->acc);
                }
                MFREE(s);
        }
}


int get_mapqual_bins(struct mapqual_bins **map)
{
        struct mapqual_bins* m = NULL;

        MMALLOC(m, sizeof(struct mapqual_bins));

        m->map = NULL;
        m->description = NULL;
        m->len = 256;
        m->n_bin = 4;
        galloc(&m->map,m->len);
        galloc(&m->description, m->n_bin, 256);

        for(int i = 0; i < m->len;i++){
                /* if(i == 0){ */
                /*         m->map[i] = MAPQUALBIN_UNMAP; */
                /* }else  */if(i < 10){
                        m->map[i] = MAPQUALBIN_lt10;
                }else if(i < 30){
                        m->map[i] = MAPQUALBIN_lt30;
                }else{
                        m->map[i] = MAPQUALBIN_gt30;
                }
        }

        snprintf(m->description[MAPQUALBIN_gt30], 256, "MAPQ gt 30");
        snprintf(m->description[MAPQUALBIN_lt30], 256, "MAPQ lt 30");
        snprintf(m->description[MAPQUALBIN_lt10], 256, "MAPQ lt 10");
        snprintf(m->description[MAPQUALBIN_UNMAP], 256, "Unmapped");


        *map = m;
        return OK;
ERROR:
        free_mapqual_bins(m);
        return FAIL;
}

void free_mapqual_bins(struct mapqual_bins *m)
{
        if(m){
                gfree(m->description);
                gfree(m->map);
                MFREE(m);
        }

}

