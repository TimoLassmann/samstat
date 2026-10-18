#ifndef REPORT_PLOTS_H
#define REPORT_PLOTS_H

/*
 * The six standard plots and the MAPQ swatches of the SAMStat report
 * (docs/plan-v2.4/standard-plots.md sections 2 and 4), drawn from a finalised
 * stat_collection on top of the svg_plot primitives.
 *
 * Output contract (the page code, WP4.2a, relies on it):
 *   - Only static markup: <svg> elements, and <p class="plot-na">,
 *     <p class="plot-note"> and <p class="plot-key"> (wraps the base key of plots 4 and 5) lines. No <script>, no <style>, no external
 *     reference, no per-point elements.
 *   - Every <svg> has role="img", <title> and <desc>. Ids are unique on the
 *     page: "samstat-p<plot>[-r1|-r2|-key]" and "samstat-mapq-<class>".
 *   - Plot 4 has one row of panels per MAPQ class with reads (aligned input) or one row of all reads;
 *     the mismatch rows of plot 5 (one per mapped class) are filled areas, one <path> per read base,
 *     stacked bottom to top A C G T N and drawn top first (standard-plots.md encoding rule 9).
 *   - Every <svg> carries style="max-width:100%;height:auto;vertical-align:top"
 *     and its natural width in the width attribute (400 units, text 15 units, so the text is
 *     still >= 12 px on a 360 px screen). A paired input is drawn as two sibling <svg>
 *     elements (Read 1, Read 2) with identical axes. Without page CSS they sit in the normal
 *     inline flow; the report page sets width:100% on them inside a grid (stat_report.c).
 *   - A plot that has no usable data is replaced by one <p class="plot-na">
 *     line (the wording lives in report_help.c). A plot that does not apply to
 *     the input type (plots 5 and 6 on unaligned input) writes nothing.
 *   - <p class="plot-note"> carries the partial-MD note of plot 5.
 *   - The output depends only on the collection: deterministic bytes.
 *
 * The caption, the "What good looks like:" line, the "How to read plot N" block and the plot title are the
 * caller's (report_help.c); these functions draw the figure only.
 */

#include "tldlite.h"

struct stat_collection;

/*
 * Palettes. The MAPQ palette is Okabe-Ito (blue, vermillion) with the orange darkened to
 * #B86E00 and a grey; every colour has at least 3:1 contrast on white (WCAG 1.4.11), the base
 * colours at least 5:1 (the thin dashed lines of plot 4 need the margin). Every
 * class and every base also has its own dash pattern, so the lines stay distinguishable
 * in greyscale and for colour-blind readers. The base palette shares no hex
 * value with the MAPQ palette (standard-plots.md 2.8). A dash of "" is solid.
 * Index = MAPQ class (COLLECT_MAPQ_*) or base (A C G T N).
 */
static const char* const report_mapq_colour[4] = {"#0072B2", "#B86E00", "#D55E00", "#777777"};
static const char* const report_mapq_dash[4] = {"", "8 4", "3 4", "12 4 2 4"};
static const char* const report_base_colour[5] = {"#007A5A", "#A8467F", "#000000", "#8C5A00", "#1F6FA8"};
static const char* const report_base_dash[5] = {"", "6 3", "2 3", "8 3 2 3", "12 3"};

/* Line of unaligned input (one line, "all reads"). */
#define REPORT_ALL_READS_COLOUR "#000000"

/*
 * Draw plot plot_id (1 = length, 2 = read quality, 3 = quality by position,
 * 4 = base composition, 5 = errors by position, 6 = identity) into o, or the
 * N/A line. Returns OK or FAIL (bad plot_id, out of memory).
 */
int report_plot_render(tld_strbuf* o, const struct stat_collection* s, int plot_id);

/*
 * The shape view of plots 1, 2 and 6 on aligned input: the same histogram with every MAPQ class
 * scaled to 100 % of its own reads (ids "samstat-p<plot>s[-r1|-r2]"). Writes nothing when it does
 * not apply: another plot, unaligned input, or fewer than two classes with 50 reads or more.
 */
int report_plot_render_shape(tld_strbuf* o, const struct stat_collection* s, int plot_id);

/* One small <svg> with the line sample (colour and dash) of a MAPQ class 0..3 (the MAPQ table). */
int report_mapq_swatch(tld_strbuf* o, int mapq_class);

#endif
