#ifndef REPORT_HELP_H
#define REPORT_HELP_H

/* Single source of all user-facing wording of the standard report
   (docs/plan-v2.4/standard-plots.md sections 1.3, 1.4, 4 and 5).
   Every string is a static literal; nothing here allocates.

   Note templates are printf format strings. The argument types are given
   next to each template; numbers are pre-formatted by the caller as strings
   so that the HTML formatting rules (thousands separators, SI suffixes) are
   applied in one place. */

#include <stddef.h>

#ifdef REPORT_HELP_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

/* Plot ids, in page order. */
enum report_plot_id {
        REPORT_PLOT_LENGTH   = 0,   /* how long are the reads? */
        REPORT_PLOT_READQ    = 1,   /* how good is each read? */
        REPORT_PLOT_POSQ     = 2,   /* where along the read does quality drop? */
        REPORT_PLOT_COMP     = 3,   /* is base composition even? */
        REPORT_PLOT_ERRORS   = 4,   /* where do errors occur? */
        REPORT_PLOT_IDENTITY = 5,   /* how accurate are the aligned reads? */
        REPORT_N_PLOTS       = 6
};

/* MAPQ classes, in table order (ids are the frozen JSON ids). */
enum report_mapq_id {
        REPORT_MAPQ_GE30     = 0,
        REPORT_MAPQ_10_29    = 1,
        REPORT_MAPQ_LT10     = 2,
        REPORT_MAPQ_UNMAPPED = 3,
        REPORT_N_MAPQ        = 4
};

struct report_plot_help {
        const char* title;          /* the question */
        const char* caption;        /* always visible */
        const char* caption_unaligned;  /* caption for input that is not aligned, or NULL for the same */
        const char* good;           /* visible line under the figure */
        const char* good_unaligned; /* good line for input that is not aligned, or NULL for the same */
        const char* look;           /* "How to read plot N" (collapsed), line 1 */
        const char* look_unaligned; /* look line for input that is not aligned, or NULL for the same */
        const char* benign;         /* line 2 */
        const char* note;           /* line 3, how a value is calculated, or NULL */
        const char* na_no_quality;  /* plots 2-3 without qualities, else NULL */
        const char* na_no_mapped;   /* plots 5-6 without mapped reads, else NULL */
        const char* na_no_data;     /* plot 6: no usable NM / =,X, else NULL */
        const char* na_no_sequence; /* plots 1, 4, 5: reads exist but none has a sequence, else NULL */
};

struct report_mapq_help {
        const char* id;             /* frozen JSON id */
        const char* label;
        const char* colour;         /* hex colour */
        const char* dash;           /* SVG stroke-dasharray, "" for solid */
};

EXTERN const struct report_plot_help report_plot_help[REPORT_N_PLOTS];
EXTERN const struct report_mapq_help report_mapq_help[REPORT_N_MAPQ];

/* Labels. */
EXTERN const char* const report_how_to_read_label;   /* "How to read plot" (the page appends the number) */
EXTERN const char* const report_note_label;          /* "How it is calculated:" */
EXTERN const char* const report_good_label;          /* "What good looks like:" */
EXTERN const char* const report_look_label;          /* "Worth a look:" */
EXTERN const char* const report_benign_label;        /* "Often benign:" */
EXTERN const char* const report_shape_label;         /* "Compare the MAPQ classes by shape" (summary of the shape view) */
EXTERN const char* const report_shape_caption;       /* what the shape view shows, and what it does not */
EXTERN const char* const report_mapq_title;          /* MAPQ table heading */
EXTERN const char* const report_mapq_caption;        /* MAPQ table caption */

/* Caption, "What good looks like:" and "Worth a look:" lines of a plot for aligned or unaligned input
   (plot 4 on unaligned input has one row of all reads, so its text names no MAPQ class). */
const char* report_plot_caption(int plot_id, int is_aligned);
const char* report_plot_good(int plot_id, int is_aligned);
const char* report_plot_look(int plot_id, int is_aligned);

/* Plot 5 mismatch row text when no read has MD or X operations. */
EXTERN const char* const report_mismatch_row_na;

/* Notes. Argument types in the comment. */
EXTERN const char* const report_note_partial;        /* (const char* N) */
EXTERN const char* const report_note_no_md;          /* (const char* N, const char* M) */
EXTERN const char* const report_note_identity_na;    /* (const char* N, const char* M) */
EXTERN const char* const report_note_beyond_axis;    /* (const char* P) */
/* Appended to the plot 5 mismatch row caption when some reads lack MD. */
EXTERN const char* const report_note_mismatch_excludes; /* (const char* N, const char* M) */

/* No sequence (SEQ is * in every read), few reads, and an empty file. */
EXTERN const char* const report_note_no_sequence;
EXTERN const char* const report_note_few_reads;      /* (const char* N, const char* MIN) */
EXTERN const char* const report_note_hidden;         /* (const char* N, const char* "s" or ""): fewer than MIN reads, plots 3-5 hidden */
EXTERN const char* const report_na_too_few;          /* (const char* N, const char* MIN): replaces plots 3-5 */
EXTERN const char* const report_banner_no_reads;
/* Sub-labels of two summary tiles. */
EXTERN const char* const report_tile_n50_sub;
EXTERN const char* const report_tile_clipped_sub;
EXTERN const char* const report_tile_meanq_sub;
EXTERN const char* const report_tile_mapq_sub;
EXTERN const char* const report_tile_identity_sub;

/* Reads below REPORT_FEW_READS get the "few reads" note. A position backed by fewer than
   REPORT_MIN_POS_BASES bases (on the absolute axis about one base per read) is a gap in plots 3-5,
   and an input with fewer reads than that gets no per-position plots at all (plots 3-5 show
   report_na_too_few): the rule that is stated is the rule that is applied. */
#define REPORT_FEW_READS 100
#define REPORT_MIN_POS_BASES 50
#define REPORT_MIN_POS_BASES_STR "50"

/* Header line suffix under -p. */
EXTERN const char* const report_header_partial;      /* (const char* N) */

#undef REPORT_HELP_IMPORT
#undef EXTERN

#endif
