#ifndef SVG_PLOT_H
#define SVG_PLOT_H

/*
 * SVG primitives for the SAMStat report (docs/plan-v2.4, WP1.6).
 *
 * Static, dependency-free inline SVG written into a tld_strbuf:
 *   - nice ticks, linear and log-10 axes
 *   - line series (dash pattern, gaps for missing values)
 *   - step histograms
 *   - one vertical marker line with a label
 *   - swatches (single, and a row of them)
 *   - panel grids (rows x columns, shared x per column, shared y per row)
 *
 * Rules (standard-plots.md section 2): every <svg> carries role="img", a
 * <title> and a <desc>; every series is ONE <path> with ONE series-level
 * <title>; there are no per-point elements, no tooltips, no direct line
 * labels and no JavaScript. The output is deterministic: the same input
 * gives the same bytes, independent of locale and thread count.
 *
 * Sizes are in SVG user units, which equal CSS px when the <svg> is shown at
 * its natural width. The root carries max-width:100%, so a figure wider than
 * its container shrinks together with its text. Callers that need >= 12 px
 * text at a 360 px viewport (328 px of content) should therefore keep
 * width * 12 / font <= 328 (the report uses 400 units with 15 unit text), or
 * stack the figures. A page may also scale a figure up with CSS width:100%.
 * A histogram bin with no drawn neighbour is at least 6 units wide and 1 unit
 * high (svg_step), so a lone bin on a log axis is not a hairline.
 *
 * Colours are passed as "#rrggbb" strings and dash patterns as
 * stroke-dasharray strings ("6 3"); NULL or "" means solid. The caller owns
 * the palette (the MAPQ and base-composition palettes live in report code).
 */

#include "tld.h"

#define SVG_DEFAULT_FONT 12.0   /* px, the minimum allowed text size */
#define SVG_MAX_TICKS 96

/* One axis. For a log axis lo must be > 0. */
typedef struct svg_axis {
        double lo;
        double hi;
        int log;                /* 0 = linear, 1 = log10 */
        const char* label;      /* axis title, may be NULL */
        int whole;              /* 1: a count, length or position - linear ticks are whole numbers */
} svg_axis;

/*
 * One data series.
 *   svg_line: n points at x[i] (or x0 + i * dx when x is NULL).
 *   svg_step: n bins; x holds n + 1 bin edges (or x0 + i * dx, i = 0..n).
 * y[i] = NaN is a gap. In svg_step a bin with y <= the axis minimum is also
 * not drawn, so empty stretches do not paint over the x axis.
 */
typedef struct svg_series {
        const char* title;      /* series-level <title>, may be NULL */
        const char* color;      /* "#0072B2" */
        const char* dash;       /* NULL or "" = solid */
        const double* x;
        double x0;
        double dx;
        const double* y;
        int n;
} svg_series;

/* Pixel box of one panel's plotting area plus its axes. Filled by svg_grid_panel. */
typedef struct svg_panel {
        double left;
        double top;
        double w;
        double h;
        double font;
        svg_axis x;
        svg_axis y;
} svg_panel;

/*
 * Rows x columns of panels in one <svg>. All plotting areas have the same
 * size. The x axis (tick labels and title) is drawn on the bottom row only
 * and the y axis (tick labels and title) on the first column only, so panels
 * must be given the same x range per column and the same y range per row.
 */
typedef struct svg_grid {
        int rows;
        int cols;
        double width;
        double height;
        double font;
        int has_titles;         /* reserve room for a title above row 0 */
        double ml;              /* margins and gaps, computed by svg_grid_init */
        double mr;
        double mt;
        double mb;
        double gap_x;
        double gap_y;
        double pw;              /* plotting area size */
        double ph;
} svg_grid;

typedef struct svg_swatch_item {
        const char* label;
        const char* color;
        const char* dash;
} svg_swatch_item;

/* ---- ticks ---- */

/* Tick step of the form 1, 2 or 5 times a power of ten giving about target ticks over range. */
double svg_nice_step(double range, int target);

/* Smallest multiple of svg_nice_step(max, target) that is >= max (max <= 0 gives 1). */
double svg_nice_ceil(double max, int target);

/*
 * Tick positions for an axis, about target labelled ticks. major[i] = 1 for a
 * labelled tick, 0 for a short unlabelled one (log axes only). Returns the
 * number of ticks (at most max), or -1 on bad input. *step is the linear tick
 * spacing (0 for log axes) and is what svg_tick_label needs.
 */
int svg_ticks(const svg_axis* a, int target, double* v, int* major, int max, double* step);

/*
 * Label for a tick: plain below 10 000 on linear axes, then "20 k" and
 * "1.5 M"; log axes use "100", "1 k", "10 k", "1 M" per tick. The label is
 * locale independent.
 */
void svg_tick_label(const svg_axis* a, double step, double v, char* buf, size_t n);

/* ---- drawing ---- */

/*
 * Open an <svg> of w x h units with base text size font (>= 12). id must be
 * unique on the page: it names the <title> and <desc> (aria-labelledby).
 */
int svg_begin(tld_strbuf* o, const char* id, double w, double h, double font,
              const char* title, const char* desc);
int svg_end(tld_strbuf* o);

void svg_grid_init(svg_grid* g, int rows, int cols, double w, double h, double font, int has_titles);

/*
 * Draw panel (r, c): grid lines, axes, ticks, tick labels and, when given, the
 * panel title (above the plot area). Fills *p for the series calls below.
 */
int svg_grid_panel(tld_strbuf* o, const svg_grid* g, int r, int c, const char* title,
                   const svg_axis* x, const svg_axis* y, svg_panel* p);

/* As svg_grid_panel; x ticks above x_tick_max are dropped (0 = no limit). */
int svg_grid_panel_ex(tld_strbuf* o, const svg_grid* g, int r, int c, const char* title,
                      const svg_axis* x, const svg_axis* y, double x_tick_max, svg_panel* p);

int svg_line(tld_strbuf* o, const svg_panel* p, const svg_series* s);
int svg_step(tld_strbuf* o, const svg_panel* p, const svg_series* s);

/* Dashed black vertical line at data value x with a text label beside it. */
int svg_vline(tld_strbuf* o, const svg_panel* p, double x, const char* label);

/* Short text centred in the plotting area of panel p (for example "none detected" in a flat panel). */
int svg_panel_note(tld_strbuf* o, const svg_panel* p, const char* text);

/* A complete small <svg> showing one line sample (for the MAPQ table). */
int svg_swatch(tld_strbuf* o, const char* id, const char* color, const char* dash,
               const char* title, const char* desc);

/* A complete <svg> with one line sample and label per item, left to right (plot 4 legend). */
int svg_swatch_row(tld_strbuf* o, const char* id, const svg_swatch_item* items, int n,
                   double font, const char* title, const char* desc);

/*
 * A key inside a figure: one line sample and label per item, left to right, wrapping onto a new
 * row when the next item would pass width. svg_key_height gives the height the key needs (rows
 * of 1.5 font); svg_key draws it with its top-left corner at (left, top).
 */
double svg_key_height(const svg_swatch_item* items, int n, double font, double width);
int svg_key(tld_strbuf* o, double left, double top, double width, const svg_swatch_item* items, int n, double font);

#endif
