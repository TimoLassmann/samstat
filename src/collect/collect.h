#ifndef COLLECT_H
#define COLLECT_H

#include <stdint.h>

#ifdef COLLECT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif
struct tl_seq_buffer;

/* MAPQ classes (m). The ids are frozen: they are the JSON ids mapq_ge30, mapq_10_29,
 * mapq_lt10 and unmapped. MAPQ 255 ("unavailable") counts as >= 30, as in 2.3. */
#define COLLECT_MAPQ_GE30 0
#define COLLECT_MAPQ_10_29 1
#define COLLECT_MAPQ_LT10 2
#define COLLECT_MAPQ_UNMAPPED 3
#define COLLECT_N_MAPQ 4

#define COLLECT_N_MATE 2           /* r: read 1 (and unpaired), read 2 */
#define COLLECT_N_BASE 5           /* A C G T N */

/* standard-plots.md 3.2: both position axes are filled for every read. */
#define COLLECT_ABS_BINS 1000      /* bin = position, positions >= 1000 are not counted */
#define COLLECT_REL_BINS 100       /* bin = (100 * position) / aligned length */

/* standard-plots.md 3.1: exact lengths 0..999, then the 1% grid g(1000) up to the first edge
 * >= 10,000,000, whose bin is the overflow bin. */
#define COLLECT_N_LEN_BINS 1922
#define COLLECT_LEN_EXACT 1000     /* lengths below this have a bin of their own */

/* standard-plots.md 3.3: per-read mean Q bins 0..93 and the error table. */
#define COLLECT_N_Q 94
#define COLLECT_QUAL_ONE 1073741824   /* 2^30: E[0], the fixed-point scale of the error table */

/* standard-plots.md 3.4: identity in 0.1% bins; bin 1000 is exactly 100%. */
#define COLLECT_N_IDENT 1001

/* Layout (3.2): "long" if the median read length is >= 1000. */
#define COLLECT_LAYOUT_SHORT 0
#define COLLECT_LAYOUT_LONG 1
#define COLLECT_LAYOUT_MIN_LONG_MEDIAN 1000

/* 128-bit unsigned integer (GCC and Clang; __extension__ keeps -pedantic quiet) */
__extension__ typedef unsigned __int128 collect_u128;

/* A value that does not apply (JSON null, left out of the report). */
#define COLLECT_NA (__builtin_nan(""))
#define COLLECT_IS_NA(x) ((x) != (x))

struct mapqual_bins {
        uint8_t* map;
        char** description;
        int len;
        int n_bin;
};

/* All accumulators of one position bin for one MAPQ class and one mate (3.2). One struct per
 * bin, so that a base touches one or two cache lines. */
struct collect_pos_bin {
        uint64_t comp[COLLECT_N_BASE];   /* read bases A C G T N; anything else counts as N */
        uint64_t qual_bases;             /* bases with a quality value */
        uint64_t qual_error_sum;         /* sum of E[q] over those bases */
        uint64_t aligned;                /* aligned read bases (M, =, X, I) of mapped reads */
        uint64_t aligned_md;             /* the same, only for reads with MD or X operations */
        uint64_t mismatches;
        uint64_t insertions;             /* insertion events (one per gap opening) */
        uint64_t deletions;              /* deletion events (one per gap opening) */
};

/* Fixed-size integer accumulators, independent of the read length (standard-plots.md 3).
 * No pointers: merging adds the arrays element by element. First index r (mate), second m
 * (MAPQ class). */
struct collect_acc {
        struct collect_pos_bin pos_abs[COLLECT_N_MATE][COLLECT_N_MAPQ][COLLECT_ABS_BINS];
        struct collect_pos_bin pos_rel[COLLECT_N_MATE][COLLECT_N_MAPQ][COLLECT_REL_BINS];
        uint64_t len_reads[COLLECT_N_MATE][COLLECT_N_MAPQ][COLLECT_N_LEN_BINS];
        uint64_t len_bases[COLLECT_N_MATE][COLLECT_N_MAPQ][COLLECT_N_LEN_BINS];
        uint64_t read_q[COLLECT_N_MATE][COLLECT_N_MAPQ][COLLECT_N_Q];
        uint64_t identity[COLLECT_N_MATE][COLLECT_N_MAPQ][COLLECT_N_IDENT];
};

struct stat_collection {
        struct mapqual_bins* mapq_map;

        /* v2 accumulators (standard-plots.md section 3) */
        struct collect_acc* acc;

        /* reads per mate and MAPQ class: basic_nums[r][m] */
        uint64_t** basic_nums;
        uint64_t n_total;          /* set by finalise: primary records */
        uint64_t n_read1;          /* set by finalise */
        uint64_t n_read2;          /* set by finalise */
        uint64_t n_paired;
        uint64_t n_proper_paired;
        uint64_t n_secondary;      /* secondary alignments (0x100), not counted as reads */
        uint64_t n_supplementary;  /* supplementary alignments (0x800), not counted as reads */

        /* scalars of 3.6 (counted per record, summed by merge) */
        uint64_t total_bases;      /* full read length, clipped bases included */
        uint64_t min_length;       /* UINT64_MAX until a read is seen; 0 after finalise if none */
        uint64_t max_length;
        uint64_t mapped_bases;     /* bases of mapped primary reads */
        uint64_t clipped_bases;    /* soft-clipped bases of mapped primary reads */
        uint64_t qual_bases_all;   /* bases with a quality value, clipped bases included */
        collect_u128 qual_error_sum_all;   /* sum of E[q] over them; internal, not in JSON */
        uint64_t reads_without_md; /* mapped reads with neither MD nor X operations */
        uint64_t identity_reads;   /* reads binned in the identity histogram */
        uint64_t identity_na;      /* mapped reads whose identity cannot be computed */
        uint64_t bases_beyond_position_axis;   /* aligned-part bases at positions >= 1000 */
        uint8_t has_sq;            /* the header has at least one @SQ line (set by the caller) */
        uint8_t is_partial_report;

        /* set by stat_collection_finalise(), 3.7; floating point values are COLLECT_NA
         * when they do not apply */
        uint8_t is_aligned;
        uint8_t layout;            /* COLLECT_LAYOUT_SHORT or COLLECT_LAYOUT_LONG */
        uint64_t n50;              /* lower edge of the storage length bin */
        uint64_t median_length;    /* lower edge of the storage length bin */
        /* The same two values for display (HTML, MultiQC, stderr line). They equal n50 and
         * median_length except where the exact scalars say better: a bin that holds one read
         * gives that read's length, and a bin edge is never shown below min_length. *_exact is
         * set when the value is the true length (below 1 kb, one read in the bin, or all reads
         * the same length), so that "~" is left off. JSON keeps the bin edges. */
        uint64_t n50_best;
        uint64_t median_best;
        uint8_t n50_exact;
        uint8_t median_exact;
        double mean_q;
        double pct_mapped;
        double pct_mapq_ge30;
        double pct_clipped;
        double median_identity;
};

/* Add the records of <sb> to <s>. Every decision is made per record. */
EXTERN int collect_stats(struct tl_seq_buffer *sb, struct stat_collection *s);
/* Add all counts of <src> to <dst> (sum; has_sq and is_partial_report OR-ed; min_length
 * minimum, max_length maximum). A collection that saw no read is the identity element.
 * Call before finalise, on collections that are not finalised. */
EXTERN int stat_collection_merge(struct stat_collection *dst, const struct stat_collection *src);
/* Derive every summary value from the accumulators (3.7). May be called more than once. */
EXTERN int stat_collection_finalise(struct stat_collection *s);
EXTERN int stat_collection_alloc(struct stat_collection **stats);
EXTERN void stat_collection_free(struct stat_collection *s);

/* Tables generated at build time (convert_tables.c). */
EXTERN uint32_t collect_qual_error(int q);                /* E[q] = round(10^(-q/10) * 2^30), q 0..93 */
EXTERN int collect_length_bin(uint64_t len);              /* storage length bin of a read length */
EXTERN uint64_t collect_length_edge(int bin);             /* lower edge of a storage length bin */
EXTERN int collect_n_display_bins(void);                  /* display length bins (plot 1 long, MultiQC) */
EXTERN uint64_t collect_display_edge(int display_bin);    /* lower edge of a display bin */
EXTERN int collect_display_bin(int storage_bin);          /* display bin that contains a storage bin */
/* Mean Q of a bin or a total from integer sums: -10 log10(error_sum / bases / 2^30);
 * COLLECT_NA if there are no bases. */
EXTERN double collect_mean_q(collect_u128 error_sum, uint64_t bases);

#undef COLLECT_IMPORT
#undef EXTERN


#endif
