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
        const char* good;           /* "How to read this", line 1 */
        const char* look;           /* line 2 */
        const char* benign;         /* line 3 */
        const char* na_no_quality;  /* plots 2-3 without qualities, else NULL */
        const char* na_no_mapped;   /* plots 5-6 without mapped reads, else NULL */
        const char* na_no_data;     /* plot 6: no usable NM / =,X, else NULL */
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
EXTERN const char* const report_how_to_read_label;   /* "How to read this" */
EXTERN const char* const report_good_label;          /* "Good:" */
EXTERN const char* const report_look_label;          /* "Worth a look:" */
EXTERN const char* const report_benign_label;        /* "Often benign:" */
EXTERN const char* const report_mapq_title;          /* MAPQ table heading */
EXTERN const char* const report_mapq_caption;        /* MAPQ table caption */

/* Plot 5 mismatch row text when no read has MD or X operations. */
EXTERN const char* const report_mismatch_row_na;

/* Notes. Argument types in the comment. */
EXTERN const char* const report_note_partial;        /* (const char* N) */
EXTERN const char* const report_note_no_md;          /* (const char* N, const char* M) */
EXTERN const char* const report_note_identity_na;    /* (const char* N, const char* M) */
EXTERN const char* const report_note_beyond_axis;    /* (const char* P) */
/* Appended to the plot 5 mismatch row caption when some reads lack MD. */
EXTERN const char* const report_note_mismatch_excludes; /* (const char* N, const char* M) */

/* Header line suffix under -p. */
EXTERN const char* const report_header_partial;      /* (const char* N) */

#undef REPORT_HELP_IMPORT
#undef EXTERN

#endif
