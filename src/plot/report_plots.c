#include "report_plots.h"

#include "svg_plot.h"
#include "../collect/collect.h"
#include "../report/report_help.h"

#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
 * Plot renderers (WP4.2b). Every plot is built the same way:
 *   1. collect the numbers of each mate into double arrays (a bump allocator, one
 *      allocation per call),
 *   2. choose the axes ONCE from all mates, so that Read 1 and Read 2 share them,
 *   3. draw each mate as its own <svg> (rows of panels sharing the x axis).
 * Draw order is unmapped, MAPQ < 10, 10-29, >= 30, so the most reliable class is on top.
 */

#define TRY(EXP) do { if ((EXP) != OK) { goto ERROR; } } while (0)

/* A single figure is FIG_W_ONE units wide and fills the page column (about 800 px, so the 14 unit
 * text is about 17 px); each figure of a Read 1 / Read 2 pair is FIG_W_PAIR units wide and fills half
 * the column. On a narrow screen a single figure keeps a minimum width and scrolls sideways
 * (stat_report.c), so that its text stays readable. */
#define FIG_W_ONE 640.0
#define FIG_W_PAIR 400.0
#define FONT 14.0
#define FONT_PAIR 15.0
#define PH_ONE 215.0      /* plotting height of a single-row figure */
#define PH_ROW 96.0       /* plotting height of one row of plot 5 */

#define MAXL 5            /* lines per panel (rule 3) */
#define MAXB 2000         /* more than the 1922 storage length bins */
#define NA_NO_READS "No reads."

/* ------------------------------------------------------------------ */
/* small helpers                                                        */
/* ------------------------------------------------------------------ */

struct pool {
        double* d;
        size_t cap;
        size_t used;
};

static double* pool_alloc(struct pool* p, size_t n)
{
        if (p->used + n > p->cap) {
                return NULL;
        }
        double* r = p->d + p->used;
        p->used += n;
        return r;   /* zero filled by calloc */
}

static int html_text(tld_strbuf* o, const char* s)
{
        const char* run = s;
        for (const char* c = s;; c++) {
                const char* rep = NULL;
                switch (*c) {
                case 0:
                        if (c > run) {
                                char* t = strndup(run, (size_t)(c - run));
                                if (!t) {
                                        return FAIL;
                                }
                                int rc = tld_append(o, t);
                                free(t);
                                return rc;
                        }
                        return OK;
                case '&':
                        rep = "&amp;";
                        break;
                case '<':
                        rep = "&lt;";
                        break;
                case '>':
                        rep = "&gt;";
                        break;
                case '"':
                        rep = "&quot;";
                        break;
                default:
                        break;
                }
                if (rep) {
                        if (c > run) {
                                char* t = strndup(run, (size_t)(c - run));
                                if (!t) {
                                        return FAIL;
                                }
                                int rc = tld_append(o, t);
                                free(t);
                                if (rc != OK) {
                                        return FAIL;
                                }
                        }
                        TRY(tld_append(o, (char*)rep));
                        run = c + 1;
                }
        }
ERROR:
        return FAIL;
}

/* <p class="..."> text </p> */
static int html_para(tld_strbuf* o, const char* cls, const char* lead, const char* text)
{
        TRY(tld_append(o, "<p class=\""));
        TRY(tld_append(o, (char*)cls));
        TRY(tld_append(o, "\">"));
        if (lead) {
                TRY(tld_append(o, "<strong>"));
                TRY(html_text(o, lead));
                TRY(tld_append(o, "</strong> "));
        }
        TRY(html_text(o, text));
        TRY(tld_append(o, "</p>\n"));
        return OK;
ERROR:
        return FAIL;
}

static int na_line(tld_strbuf* o, const char* text)
{
        return html_para(o, "plot-na", NULL, text);
}

/* 1234567 -> "1,234,567" */
static void fmt_count(char* buf, size_t n, uint64_t v)
{
        char raw[32];
        char out[48];
        int len = snprintf(raw, sizeof(raw), "%" PRIu64, v);
        int k = 0;
        for (int i = 0; i < len; i++) {
                if (i > 0 && (len - i) % 3 == 0) {
                        out[k++] = ',';
                }
                out[k++] = raw[i];
        }
        out[k] = 0;
        snprintf(buf, n, "%s", out);
}

/* Plots 3-5 of an input with fewer reads than REPORT_MIN_POS_BASES: single reads, not a profile. */
static int too_few_line(tld_strbuf* o, const struct stat_collection* s)
{
        char n[48], text[200];
        fmt_count(n, sizeof(n), s->n_total);
        snprintf(text, sizeof(text), report_na_too_few, n, REPORT_MIN_POS_BASES_STR);
        return na_line(o, text);
}

/* Locale independent "%.<d>f". */
static void fmt_dec(char* buf, size_t n, int decimals, double v)
{
        snprintf(buf, n, "%.*f", decimals, v);
        for (char* c = buf; *c; c++) {
                if (*c == ',') {
                        *c = '.';
                }
        }
}

/* Lengths as in the summary strip: bp below 10 kb, kb up to 1 Mb, Mb above. A value above 1 kb that is a bin
 * edge (within 1%) is marked with "~"; exact is set when the value is the true length. */
static void fmt_len(char* buf, size_t n, uint64_t v, int exact)
{
        char t[48];
        const char* tilde = (v < 1000 || exact) ? "" : "~";
        if (v < 1000) {
                fmt_count(t, sizeof(t), v);
                snprintf(buf, n, "%s bp", t);
        } else if (v < 10000) {
                fmt_count(t, sizeof(t), v);
                snprintf(buf, n, "%s%s bp", tilde, t);
        } else if (v >= 999950) {
                char d[32];
                fmt_dec(d, sizeof(d), 1, (double)v / 1000000.0);
                snprintf(buf, n, "%s%s Mb", tilde, d);
        } else {
                char d[32];
                fmt_dec(d, sizeof(d), 1, (double)v / 1000.0);
                /* thousands separators in the integer part of the kb value */
                char* dot = strchr(d, '.');
                size_t il = dot ? (size_t)(dot - d) : strlen(d);
                char ip[32];
                snprintf(ip, sizeof(ip), "%.*s", (int)il, d);
                fmt_count(t, sizeof(t), (uint64_t)strtoull(ip, NULL, 10));
                snprintf(buf, n, "%s%s%s kb", tilde, t, dot ? dot : "");
        }
}

static int mates(const struct stat_collection* s)
{
        return s->n_read2 > 0 ? 2 : 1;
}

static uint64_t n_mapped(const struct stat_collection* s)
{
        uint64_t n = 0;
        for (int r = 0; r < COLLECT_N_MATE; r++) {
                for (int m = 0; m < COLLECT_MAPQ_UNMAPPED; m++) {
                        n += s->basic_nums[r][m];
                }
        }
        return n;
}

/* ------------------------------------------------------------------ */
/* lines                                                                */
/* ------------------------------------------------------------------ */

struct lineinfo {
        const char* label;
        const char* colour;
        const char* dash;
        unsigned mask;          /* MAPQ classes summed into this line */
};

/*
 * Lines of a MAPQ-split plot, in draw order (unmapped first, so that MAPQ >= 30 ends up on
 * top). Aligned input: one line per class (the unmapped class only with_unmapped). Unaligned
 * input: one solid black line, "all reads", all classes pooled (F7).
 */
static int get_lines(const struct stat_collection* s, int with_unmapped, struct lineinfo* l)
{
        int n = 0;
        if (!s->is_aligned) {
                l[0] = (struct lineinfo){"All reads", REPORT_ALL_READS_COLOUR, "", 0xFu};
                return 1;
        }
        for (int m = COLLECT_MAPQ_UNMAPPED; m >= 0; m--) {
                if (m == COLLECT_MAPQ_UNMAPPED && !with_unmapped) {
                        continue;
                }
                l[n++] = (struct lineinfo){report_mapq_help[m].label, report_mapq_colour[m], report_mapq_dash[m], 1u << m};
        }
        return n;
}

/* ------------------------------------------------------------------ */
/* generic figure: rows of panels for one mate                          */
/* ------------------------------------------------------------------ */

struct panel_data {
        const char* title;
        svg_axis x;
        svg_axis y;
        double x_tick_max;      /* x ticks above this are dropped; 0 = no limit */
        int step;               /* 1 = step histogram, 0 = line */
        int n;
        const double* xe;       /* x values / bin edges, or NULL for x0 + i * dx */
        double x0;
        double dx;
        int nl;
        struct lineinfo li[MAXL];
        const double* v[MAXL];  /* NULL = no data for this line */
        int has_vline;
        double vline;
        char vlabel[64];
        const char* note;       /* short text inside the panel (plot 5: "none detected"), or NULL */
};

/* Top of a count axis: a little headroom over the tallest bar, and never below 3 (4 from 3 reads on)
 * so that a handful of reads gets whole-number ticks, not ticks at 0.2 reads. */
static double count_top(double ymax)
{
        double t = svg_nice_ceil(ymax * 1.04, 5);
        if (t < 4.0) {
                return ymax <= 2.0 ? 3.0 : 4.0;
        }
        return t;
}

static int has_values(const double* v, int n, int step)
{
        if (!v) {
                return 0;
        }
        for (int i = 0; i < n; i++) {
                if (isfinite(v[i]) && (!step || v[i] > 0.0)) {
                        return 1;
                }
        }
        return 0;
}

/* "MAPQ >= 30, MAPQ 10-29" of the lines that have data in this figure's first row. */
static void list_lines(const struct panel_data* p, char* buf, size_t n)
{
        buf[0] = 0;
        for (int k = p->nl - 1; k >= 0; k--) {      /* table order: MAPQ >= 30 first */
                if (has_values(p->v[k], p->n, p->step)) {
                        size_t l = strlen(buf);
                        snprintf(buf + l, n - l, "%s%s", l ? ", " : "", p->li[k].label);
                }
        }
        if (!buf[0]) {
                snprintf(buf, n, "none");
        }
}

/*
 * Key of a MAPQ-split figure, in table order (MAPQ >= 30 first): every class with reads in this
 * mate, with its read count, so that a flat line at zero explains itself. A class that has reads
 * but no line in any row (plots 3 and 5: fewer than 50 reads at every position; the shape view:
 * fewer than 50 reads) says so. Returns the number of items.
 */
#define KEY_LABEL 96
static int mapq_key(const struct stat_collection* s, int mate, const struct panel_data* rows, int nrows,
                    svg_swatch_item* items, char (*labels)[KEY_LABEL])
{
        const struct panel_data* d = &rows[0];
        int n = 0;

        if (!s || !s->is_aligned || d->nl < 2) {
                return 0;
        }
        for (int k = d->nl - 1; k >= 0; k--) {
                uint64_t reads = 0, any = 0;
                char cnt[48];
                int drawn = 0;
                for (int m = 0; m < COLLECT_N_MAPQ; m++) {
                        if (d->li[k].mask & (1u << m)) {
                                reads += s->basic_nums[mate][m];
                                for (int r = 0; r < mates(s); r++) {
                                        any += s->basic_nums[r][m];
                                }
                        }
                }
                /* the classes of either mate: Read 1 and Read 2 get keys of the same height, so that
                 * their panels stay side by side */
                if (any == 0) {
                        continue;
                }
                for (int r = 0; r < nrows; r++) {
                        drawn |= has_values(rows[r].v[k], rows[r].n, rows[r].step);
                }
                fmt_count(cnt, sizeof(cnt), reads);
                snprintf(labels[n], KEY_LABEL, "%s (%s read%s%s)", d->li[k].label, cnt, reads == 1 ? "" : "s",
                         drawn || reads == 0 ? "" : ", too few to draw");
                items[n] = (svg_swatch_item){labels[n], d->li[k].colour, d->li[k].dash};
                n++;
        }
        return n;
}

/* One figure of one mate. key_s: the collection whose MAPQ classes the lines are (a key is drawn
 * above the panels), or NULL for a figure whose lines are not MAPQ classes (plot 4). */
static int figure(tld_strbuf* o, const struct stat_collection* key_s, const char* id_base, const char* title,
                  const char* desc, int mate, int nm, const struct panel_data* rows, int nrows, double ph)
{
        svg_grid g;
        svg_panel p;
        char id[64];
        char ttl[256];
        char sertitle[128];
        const int titles = 1;
        svg_swatch_item key[MAXL];
        char key_label[MAXL][KEY_LABEL];

        const double fw = nm == 2 ? FIG_W_PAIR : FIG_W_ONE;
        const double font = nm == 2 ? FONT_PAIR : FONT;
        const int nkey = mapq_key(key_s, mate, rows, nrows, key, key_label);

        svg_grid_init(&g, nrows, 1, fw, 300.0, font, titles);
        const double kh = nkey ? svg_key_height(key, nkey, font, fw - 1.5 * font) + 0.5 * font : 0.0;
        double h = kh + g.mt + g.mb + (nrows - 1) * g.gap_y + nrows * ph;
        svg_grid_init(&g, nrows, 1, fw, h, font, titles);
        /* the key sits above the first panel: move the panels down by its height */
        g.mt += kh;
        g.ph = ph;

        if (nm == 2) {
                snprintf(id, sizeof(id), "%s-r%d", id_base, mate + 1);
                snprintf(ttl, sizeof(ttl), "%s Read %d", title, mate + 1);
        } else {
                snprintf(id, sizeof(id), "%s", id_base);
                snprintf(ttl, sizeof(ttl), "%s", title);
        }
        TRY(svg_begin(o, id, fw, h, font, ttl, desc));
        if (nkey) {
                TRY(svg_key(o, 0.5 * font, 0.5 * font, fw - 1.5 * font, key, nkey, font));
        }
        for (int r = 0; r < nrows; r++) {
                const struct panel_data* d = &rows[r];
                TRY(svg_grid_panel_ex(o, &g, r, 0, d->title, &d->x, &d->y, d->x_tick_max, &p));
                for (int k = 0; k < d->nl; k++) {
                        if (!d->v[k]) {
                                continue;
                        }
                        if (nm == 2) {
                                snprintf(sertitle, sizeof(sertitle), "%s, Read %d", d->li[k].label, mate + 1);
                        } else {
                                snprintf(sertitle, sizeof(sertitle), "%s", d->li[k].label);
                        }
                        svg_series ser = {sertitle, d->li[k].colour, d->li[k].dash, d->xe, d->x0, d->dx, d->v[k], d->n};
                        TRY(d->step ? svg_step(o, &p, &ser) : svg_line(o, &p, &ser));
                }
                if (d->has_vline) {
                        TRY(svg_vline(o, &p, d->vline, d->vlabel));
                }
                if (d->note) {
                        TRY(svg_panel_note(o, &p, d->note));
                }
        }
        return svg_end(o);
ERROR:
        return FAIL;
}

/* ------------------------------------------------------------------ */
/* axes                                                                 */
/* ------------------------------------------------------------------ */

static double nice125_down(double v)
{
        double p = pow(10.0, floor(log10(v) + 1e-12));
        double m = v / p;
        return (m >= 5.0 - 1e-9 ? 5.0 : m >= 2.0 - 1e-9 ? 2.0 : 1.0) * p;
}

static double nice125_up(double v)
{
        double p = pow(10.0, floor(log10(v) + 1e-12));
        double m = v / p;
        return (m <= 1.0 + 1e-9 ? 1.0 : m <= 2.0 + 1e-9 ? 2.0 : m <= 5.0 + 1e-9 ? 5.0 : 10.0) * p;
}

struct posaxis {
        int rel;                /* 1 = relative axis (100 bins), 0 = absolute (1000 bins) */
        int nb;                 /* bins drawn */
        double x0;              /* x of the first bin */
        svg_axis x;             /* x.label points into lab_abs/lab_rel: per call, never shared between threads */
        char lab_abs[64];
        char lab_rel[64];
        uint64_t min_bases;     /* a bin backed by fewer bases than this (all classes, one mate) is a gap */
};

static uint64_t pos_bin_bases(const struct stat_collection* s, const struct posaxis* a, int r, int i)
{
        uint64_t t = 0;
        for (int m = 0; m < COLLECT_N_MAPQ; m++) {
                const struct collect_pos_bin* b = a->rel ? &s->acc->pos_rel[r][m][i] : &s->acc->pos_abs[r][m][i];
                for (int c = 0; c < COLLECT_N_BASE; c++) {
                        t += b->comp[c];
                }
        }
        return t;
}

/*
 * Reads of the classes in <mask> that cover position bin i of mate r. Absolute axis: a read puts at
 * most one base into a position, so the bases in the bin are the reads. Relative axis: every read
 * spans all 100 bins, so it is the number of reads in the classes. A line of plots 3 and 5 is drawn
 * only where its own classes have at least a->min_bases reads: a class of a handful of reads is noise,
 * and would also set the vertical scale for the classes that are measured.
 */
static uint64_t line_support(const struct stat_collection* s, const struct posaxis* a, int r, unsigned mask, int i)
{
        uint64_t t = 0;
        for (int m = 0; m < COLLECT_N_MAPQ; m++) {
                if (!(mask & (1u << m))) {
                        continue;
                }
                if (a->rel) {
                        t += s->basic_nums[r][m];
                } else {
                        const struct collect_pos_bin* b = &s->acc->pos_abs[r][m][i];
                        for (int c = 0; c < COLLECT_N_BASE; c++) {
                                t += b->comp[c];
                        }
                }
        }
        return t;
}

/*
 * x axis of plots 3-5 (3.2). Short layout: absolute, from position 1 to the last bin with any
 * base (over all classes and both mates). Long layout: relative, 0-100 %. Returns 0 when the input
 * has fewer than REPORT_MIN_POS_BASES reads or no bin is backed by enough bases: the caller shows
 * the "too few reads" line, so that the rule the page states is the rule that is applied.
 */
static int pos_axis(const struct stat_collection* s, struct posaxis* a)
{
        const int nb_max = (s->layout == COLLECT_LAYOUT_LONG) ? COLLECT_REL_BINS : COLLECT_ABS_BINS;

        a->rel = (s->layout == COLLECT_LAYOUT_LONG);
        snprintf(a->lab_abs, sizeof(a->lab_abs), "position in read (bp)%s", s->is_aligned ? ", aligned part" : "");
        snprintf(a->lab_rel, sizeof(a->lab_rel), "position in read (%% of length)%s", s->is_aligned ? ", aligned part" : "");

        /* last bin that any mate backs with enough bases */
        int last = -1;
        a->min_bases = REPORT_MIN_POS_BASES;
        a->nb = 0;
        a->x0 = 1.0;
        a->x = (svg_axis){0.5, 1.5, 0, a->lab_abs, 1};
        if (s->n_total < REPORT_MIN_POS_BASES) {
                return 0;
        }
        for (int r = 0; r < COLLECT_N_MATE; r++) {
                for (int i = nb_max - 1; i > last; i--) {
                        if (pos_bin_bases(s, a, r, i) >= a->min_bases) {
                                last = i;
                                break;
                        }
                }
        }
        if (a->rel) {
                a->nb = COLLECT_REL_BINS;
                a->x0 = 0.5;
                a->x = (svg_axis){0.0, 100.0, 0, a->lab_rel, 1};
                return last >= 0;
        }
        a->nb = last + 1;
        a->x0 = 1.0;
        /* positions 1..nb are drawn at x = 1..nb; half a bin of room on each side keeps the lines off the
         * frame and lets the ticks fall on 50, 100, 150 */
        a->x = (svg_axis){0.5, (double)(a->nb > 0 ? a->nb : 1) + 0.5, 0, a->lab_abs, 1};
        return a->nb > 0;
}

static const struct collect_pos_bin* pos_bins(const struct stat_collection* s, const struct posaxis* a, int r, int m)
{
        return a->rel ? s->acc->pos_rel[r][m] : s->acc->pos_abs[r][m];
}

/* ------------------------------------------------------------------ */
/* plot 1: how long are the reads?                                      */
/* ------------------------------------------------------------------ */

/*
 * Long layout: the 1% storage bins are far too fine to draw (about 230 per decade, so a few
 * thousand reads leave most bins empty and the outline becomes a comb). The figure merges them
 * into LEN_PER_DECADE bins per decade, each a union of whole storage bins, so the counts stay
 * exact: storage bin b goes into plot bin floor(LEN_PER_DECADE * log10(lower edge of b)).
 * Short layout: one bin per length, drawn from length - 0.5 to length + 0.5 so that the bar is
 * centred on its length.
 */
#define LEN_PER_DECADE 20

static int len_plot_bin(int is_long, int b)
{
        if (!is_long) {
                return b;
        }
        uint64_t e = collect_length_edge(b);
        return e <= 1 ? 0 : (int)floor(LEN_PER_DECADE * log10((double)e) + 1e-9);
}

static int len_n_plot_bins(int is_long)
{
        return is_long ? len_plot_bin(1, COLLECT_N_LEN_BINS - 1) + 1 : COLLECT_N_LEN_BINS;
}

/* Lower edges of the plot bins, n + 1 values: e[d] is the smallest storage edge in plot bin d or
 * above (a plot bin without a storage bin, possible below 10 bp, has zero width); e[n] is the end
 * of the overflow bin, one more 1% step. */
static void len_plot_edges(int is_long, double* e, int n)
{
        const double top = (double)collect_length_edge(COLLECT_N_LEN_BINS - 1) * 1.01 + 1.0;
        const double shift = is_long ? 0.0 : 0.5;

        for (int d = 0; d <= n; d++) {
                e[d] = top;
        }
        for (int b = COLLECT_N_LEN_BINS - 1; b >= 0; b--) {
                e[len_plot_bin(is_long, b)] = (double)collect_length_edge(b);
        }
        for (int d = n - 1; d >= 0; d--) {
                if (e[d] > e[d + 1]) {
                        e[d] = e[d + 1];
                }
        }
        for (int d = 0; d <= n; d++) {
                e[d] -= shift;
        }
}


/*
 * The shape view of plots 1, 2 and 6 (the second, collapsed figure): every line is scaled to the
 * reads of its own class, so that it adds up to 100 % and a class of 40 reads can be compared with
 * one of 40 000 by shape. A class with fewer than REPORT_MIN_POS_BASES reads in a mate is left out
 * (its array is set to NULL): a handful of reads would draw a jagged line as tall as the others.
 * cnt[r][k][off .. off + n - 1] are the drawn bins. Sets *ymax to the tallest share and returns 1
 * when some mate keeps at least two lines (one line has nothing to be compared with), else 0.
 */
static int shape_scale(double* cnt[2][MAXL], int nm, int nl, int off, int n, double* ymax)
{
        int any = 0;

        *ymax = 0.0;
        for (int r = 0; r < nm; r++) {
                int kept = 0;
                for (int k = 0; k < nl; k++) {
                        double t = 0.0;
                        if (!cnt[r][k]) {
                                continue;
                        }
                        for (int d = off; d < off + n; d++) {
                                t += cnt[r][k][d];
                        }
                        if (t < (double)REPORT_MIN_POS_BASES) {
                                cnt[r][k] = NULL;
                                continue;
                        }
                        for (int d = off; d < off + n; d++) {
                                cnt[r][k][d] *= 100.0 / t;
                                *ymax = cnt[r][k][d] > *ymax ? cnt[r][k][d] : *ymax;
                        }
                        kept++;
                }
                any |= (kept >= 2);
        }
        return any;
}

static svg_axis share_axis(double ymax)
{
        return (svg_axis){0.0, svg_nice_ceil(ymax * 1.04, 5), 0, "% of the reads in the class", 0};
}

static void shape_desc(char* buf, size_t n, const char* what, const char* names)
{
        snprintf(buf, n,
                 "Share of the reads of each MAPQ class %s. Each line adds up to 100%%, so that classes of "
                 "different sizes are compared by shape, not by number. Classes with fewer than %s reads are "
                 "left out. Lines: %s.", what, REPORT_MIN_POS_BASES_STR, names);
}

/* Index into report_plot_help[] for plot ids 1..6. */
static const struct report_plot_help* help_of(int plot_id)
{
        return &report_plot_help[plot_id - 1];
}

static const char* mate_title(int nm, int r)
{
        return nm == 2 ? (r ? "Read 2" : "Read 1") : NULL;
}

/* ------------------------------------------------------------------ */
/* plot 1: how long are the reads?                                      */
/* ------------------------------------------------------------------ */

static int plot_length(tld_strbuf* o, const struct stat_collection* s, struct pool* pl, int shape)
{
        const int is_long = (s->layout == COLLECT_LAYOUT_LONG);
        const int nb = len_n_plot_bins(is_long);
        const int nm = mates(s);
        struct lineinfo li[MAXL];
        const int nl = get_lines(s, 1, li);
        double* cnt[2][MAXL];
        int bmin = nb, bmax = -1;
        double ymax = 0.0;

        if (shape && (!s->is_aligned || s->total_bases == 0 || s->min_length == s->max_length)) {
                return OK;      /* no shape view: one line, or no histogram */
        }
        if (s->total_bases == 0) {
                return na_line(o, help_of(1)->na_no_sequence);
        }
        if (s->min_length == s->max_length) {
                /* a histogram of one length is a single bar: the sentence says it better */
                char nr[48], len[48], text[160];
                fmt_count(nr, sizeof(nr), s->n_total);
                fmt_len(len, sizeof(len), s->max_length, 1);
                if (s->n_total == 1) {
                        snprintf(text, sizeof(text), "The read is %s long.", len);
                } else {
                        snprintf(text, sizeof(text), "All %s reads are %s long.", nr, len);
                }
                return html_para(o, "plot-note", NULL, text);
        }
        for (int r = 0; r < nm; r++) {
                for (int k = 0; k < nl; k++) {
                        cnt[r][k] = pool_alloc(pl, (size_t)nb);
                        if (!cnt[r][k]) {
                                return FAIL;
                        }
                        for (int m = 0; m < COLLECT_N_MAPQ; m++) {
                                if (!(li[k].mask & (1u << m))) {
                                        continue;
                                }
                                for (int b = 0; b < COLLECT_N_LEN_BINS; b++) {
                                        uint64_t v = s->acc->len_reads[r][m][b];
                                        if (v) {
                                                cnt[r][k][len_plot_bin(is_long, b)] += (double)v;
                                        }
                                }
                        }
                        for (int d = 0; d < nb; d++) {
                                if (cnt[r][k][d] > 0.0) {
                                        bmin = d < bmin ? d : bmin;
                                        bmax = d > bmax ? d : bmax;
                                        ymax = cnt[r][k][d] > ymax ? cnt[r][k][d] : ymax;
                                }
                        }
                }
        }
        if (bmax < 0) {
                return shape ? OK : na_line(o, NA_NO_READS);
        }

        const int n = bmax - bmin + 1;
        if (shape && !shape_scale(cnt, nm, nl, bmin, n, &ymax)) {
                return OK;
        }
        double* edges = pool_alloc(pl, (size_t)nb + 1);
        if (!edges) {
                return FAIL;
        }
        len_plot_edges(is_long, edges, nb);
        double* xe = edges + bmin;      /* n + 1 edges */

        svg_axis xa;
        if (is_long) {
                double lo = fmax(xe[0], 1.0);
                double hi = xe[n];
                if (hi < lo * 10.0) {   /* at least one decade, centred on the data */
                        double mid = sqrt(lo * hi);
                        lo = mid / sqrt(10.0);
                        hi = mid * sqrt(10.0);
                }
                xa = (svg_axis){nice125_down(fmax(lo, 1.0)), nice125_up(hi), 1, "read length (bp), log scale", 0};
        } else {
                double lo = xe[0];
                double hi = xe[n];
                if (hi - lo < 10.0) {   /* at least 10 bp, centred on the data */
                        double mid = (lo + hi) / 2.0;
                        lo = mid - 5.0;
                        hi = mid + 5.0;
                        if (lo < 0.0) {
                                lo = 0.0;
                                hi = 10.0;
                        }
                }
                xa = (svg_axis){lo, hi, 0, "read length (bp)", 1};
        }
        svg_axis ya = shape ? share_axis(ymax) : (svg_axis){0.0, count_top(ymax), 0, "number of reads", 1};

        char lo_s[48], hi_s[48], n50_s[48], names[200], desc[900];
        fmt_len(lo_s, sizeof(lo_s), s->min_length, 1);
        fmt_len(hi_s, sizeof(hi_s), s->max_length, 1);
        fmt_len(n50_s, sizeof(n50_s), s->n50_best, s->n50_exact);

        for (int r = 0; r < nm; r++) {
                struct panel_data pd;
                memset(&pd, 0, sizeof(pd));
                pd.title = mate_title(nm, r);
                pd.x = xa;
                pd.y = ya;
                pd.step = 1;
                pd.n = n;
                pd.xe = xe;
                pd.nl = nl;
                for (int k = 0; k < nl; k++) {
                        pd.li[k] = li[k];
                        pd.v[k] = cnt[r][k] && has_values(cnt[r][k] + bmin, n, 1) ? cnt[r][k] + bmin : NULL;
                }
                if (!shape && s->n50_best > 0) {
                        pd.has_vline = 1;
                        pd.vline = (double)s->n50_best;
                        snprintf(pd.vlabel, sizeof(pd.vlabel), "N50 %s", n50_s);
                }
                list_lines(&pd, names, sizeof(names));
                snprintf(desc, sizeof(desc),
                         "Number of reads per read-length bin%s. Read lengths range from %s to %s. Lines: %s. "
                         "The dashed vertical line marks N50, %s: half of all bases are in reads at least this long.",
                         is_long ? " on a log scale" : "", lo_s, hi_s, names, n50_s);
                if (shape) {
                        shape_desc(desc, sizeof(desc), is_long ? "per read-length bin, on a log scale" : "per read length", names);
                }
                if (figure(o, s, shape ? "samstat-p1s" : "samstat-p1", help_of(1)->title, desc, r, nm, &pd, 1, PH_ONE) != OK) {
                        return FAIL;
                }
        }
        return OK;
}

/* ------------------------------------------------------------------ */
/* plot 2: how good is each read?                                       */
/* ------------------------------------------------------------------ */

static int plot_readq(tld_strbuf* o, const struct stat_collection* s, struct pool* pl, int shape)
{
        const int nm = mates(s);
        struct lineinfo li[MAXL];
        const int nl = get_lines(s, 1, li);
        double* cnt[2][MAXL];
        int qmax = -1;
        int qmin = COLLECT_N_Q;
        double ymax = 0.0;

        if (shape && (!s->is_aligned || s->qual_bases_all == 0)) {
                return OK;
        }
        if (s->qual_bases_all == 0) {
                return na_line(o, help_of(2)->na_no_quality);
        }
        for (int r = 0; r < nm; r++) {
                for (int k = 0; k < nl; k++) {
                        cnt[r][k] = pool_alloc(pl, COLLECT_N_Q);
                        if (!cnt[r][k]) {
                                return FAIL;
                        }
                        for (int m = 0; m < COLLECT_N_MAPQ; m++) {
                                if (!(li[k].mask & (1u << m))) {
                                        continue;
                                }
                                for (int q = 0; q < COLLECT_N_Q; q++) {
                                        cnt[r][k][q] += (double)s->acc->read_q[r][m][q];
                                }
                        }
                        for (int q = 0; q < COLLECT_N_Q; q++) {
                                if (cnt[r][k][q] > 0.0) {
                                        qmax = q > qmax ? q : qmax;
                                        qmin = q < qmin ? q : qmin;
                                        ymax = cnt[r][k][q] > ymax ? cnt[r][k][q] : ymax;
                                }
                        }
                }
        }
        if (qmax < 0) {
                return shape ? OK : na_line(o, help_of(2)->na_no_quality);
        }
        /*
         * The x axis spans the data, from the multiple of 5 below the first bin to the multiple of 5
         * past the last one (a bar is never the frame). It does not start at Q0: the counts are the
         * bar heights, which start at 0; the quality is a position, and reads between Q25 and Q35
         * squeezed into the right quarter of a Q0-Q35 axis would be hard to read.
         */
        const int n = qmax - qmin + 1;
        if (shape && !shape_scale(cnt, nm, nl, qmin, n, &ymax)) {
                return OK;
        }
        svg_axis xa = {5.0 * floor((double)(qmin - 1) / 5.0), 5.0 * floor((double)(qmax + 1) / 5.0 + 1.0), 0,
                       "average base quality of the read (Q)", 1};
        if (xa.lo < 0.0) {
                xa.lo = 0.0;
        }
        svg_axis ya = shape ? share_axis(ymax) : (svg_axis){0.0, count_top(ymax), 0, "number of reads", 1};
        char names[200], desc[700];

        for (int r = 0; r < nm; r++) {
                struct panel_data pd;
                memset(&pd, 0, sizeof(pd));
                pd.title = mate_title(nm, r);
                pd.x = xa;
                pd.y = ya;
                pd.step = 1;
                pd.n = n;
                pd.x0 = (double)qmin;
                pd.dx = 1.0;
                pd.nl = nl;
                for (int k = 0; k < nl; k++) {
                        pd.li[k] = li[k];
                        pd.v[k] = cnt[r][k] && has_values(cnt[r][k] + qmin, n, 1) ? cnt[r][k] + qmin : NULL;
                }
                list_lines(&pd, names, sizeof(names));
                snprintf(desc, sizeof(desc),
                         "Number of reads by their mean base quality, from Q%d to Q%d. Lines: %s. "
                         "Qualities are averaged as error probabilities.", qmin, qmax, names);
                if (shape) {
                        shape_desc(desc, sizeof(desc), "by mean base quality of the read", names);
                }
                if (figure(o, s, shape ? "samstat-p2s" : "samstat-p2", help_of(2)->title, desc, r, nm, &pd, 1, PH_ONE) != OK) {
                        return FAIL;
                }
        }
        return OK;
}

/* ------------------------------------------------------------------ */
/* plot 3: where along the read does quality drop?                      */
/* ------------------------------------------------------------------ */

static void pos_range_text(const struct posaxis* a, char* buf, size_t n)
{
        if (a->rel) {
                snprintf(buf, n, "from 0 to 100 percent of the read length");
        } else {
                snprintf(buf, n, "from position 1 to position %d", a->nb > 1 ? a->nb : 2);
        }
}

static int plot_posq(tld_strbuf* o, const struct stat_collection* s, struct pool* pl)
{
        const int nm = mates(s);
        struct lineinfo li[MAXL];
        const int nl = get_lines(s, 1, li);
        struct posaxis ax;
        double* v[2][MAXL];
        double qmax = 0.0;

        if (s->qual_bases_all == 0) {
                return na_line(o, help_of(3)->na_no_quality);
        }
        if (!pos_axis(s, &ax)) {
                return too_few_line(o, s);
        }
        for (int r = 0; r < nm; r++) {
                for (int k = 0; k < nl; k++) {
                        v[r][k] = pool_alloc(pl, (size_t)ax.nb);
                        if (!v[r][k]) {
                                return FAIL;
                        }
                        for (int i = 0; i < ax.nb; i++) {
                                uint64_t bases = 0;
                                collect_u128 err = 0;
                                for (int m = 0; m < COLLECT_N_MAPQ; m++) {
                                        if (li[k].mask & (1u << m)) {
                                                const struct collect_pos_bin* pb = &pos_bins(s, &ax, r, m)[i];
                                                bases += pb->qual_bases;
                                                err += (collect_u128)pb->qual_error_sum;
                                        }
                                }
                                double q = line_support(s, &ax, r, li[k].mask, i) >= ax.min_bases
                                                   ? collect_mean_q(err, bases) : NAN;   /* NaN = gap */
                                v[r][k][i] = q;
                                if (isfinite(q) && q > qmax) {
                                        qmax = q;
                                }
                        }
                }
        }
        svg_axis ya = {0.0, fmax(10.0, ceil(qmax / 10.0 - 1e-9) * 10.0), 0, "mean base quality (Q)", 1};
        char names[200], desc[800], range[96];
        pos_range_text(&ax, range, sizeof(range));

        for (int r = 0; r < nm; r++) {
                struct panel_data pd;
                memset(&pd, 0, sizeof(pd));
                pd.title = mate_title(nm, r);
                pd.x = ax.x;
                pd.y = ya;
                pd.n = ax.nb;
                pd.x0 = ax.x0;
                pd.dx = 1.0;
                pd.nl = nl;
                for (int k = 0; k < nl; k++) {
                        pd.li[k] = li[k];
                        pd.v[k] = has_values(v[r][k], ax.nb, 0) ? v[r][k] : NULL;
                }
                list_lines(&pd, names, sizeof(names));
                snprintf(desc, sizeof(desc),
                         "Mean base quality at each position along the read, %s%s. Lines: %s. "
                         "Gaps mean that no base with a quality was seen or that fewer than " REPORT_MIN_POS_BASES_STR
                         " reads of that line cover the position.",
                         range, s->is_aligned ? ", aligned part of the read" : "", names);
                if (figure(o, s, "samstat-p3", help_of(3)->title, desc, r, nm, &pd, 1, PH_ONE) != OK) {
                        return FAIL;
                }
        }
        return OK;
}

/* ------------------------------------------------------------------ */
/* plot 4: is base composition even along the read?                     */
/* ------------------------------------------------------------------ */

static const char* const base_name[5] = {"A", "C", "G", "T", "N"};
static const char* const base_series[5] = {"Base A", "Base C", "Base G", "Base T", "Base N"};
static const char* const base_colour_name[5] = {"green", "reddish purple", "black", "brown", "blue"};
static const char* const base_dash_name[5] = {"solid line", "dashed line", "dotted line", "dash-dot line", "long-dashed line"};

static int plot_comp(tld_strbuf* o, const struct stat_collection* s, struct pool* pl)
{
        const int nm = mates(s);
        struct posaxis ax;
        double* v[2][5];

        if (s->total_bases == 0) {
                return na_line(o, help_of(4)->na_no_sequence);
        }
        if (!pos_axis(s, &ax)) {
                return too_few_line(o, s);
        }
        /* the key: one row of five swatches above the panels */
        {
                svg_swatch_item items[5];
                char desc[300];
                int l = snprintf(desc, sizeof(desc), "Line key for the bases:");
                for (int b = 0; b < 5; b++) {
                        items[b] = (svg_swatch_item){base_name[b], report_base_colour[b], report_base_dash[b]};
                        l += snprintf(desc + l, sizeof(desc) - (size_t)l, "%s %s %s, %s", b ? ";" : "", base_name[b],
                                      base_colour_name[b], base_dash_name[b]);
                }
                if (tld_append(o, "<p class=\"plot-key\">") != OK ||
                    svg_swatch_row(o, "samstat-p4-key", items, 5, FONT, "Line key for the base composition plot", desc) != OK ||
                    tld_append(o, "</p>\n") != OK) {
                        return FAIL;
                }
        }
        for (int r = 0; r < nm; r++) {
                for (int b = 0; b < 5; b++) {
                        v[r][b] = pool_alloc(pl, (size_t)ax.nb);
                        if (!v[r][b]) {
                                return FAIL;
                        }
                }
                for (int i = 0; i < ax.nb; i++) {
                        uint64_t c[5] = {0, 0, 0, 0, 0}, tot = 0;
                        for (int m = 0; m < COLLECT_N_MAPQ; m++) {
                                const struct collect_pos_bin* pb = &pos_bins(s, &ax, r, m)[i];
                                for (int b = 0; b < 5; b++) {
                                        c[b] += pb->comp[b];
                                }
                        }
                        for (int b = 0; b < 5; b++) {
                                tot += c[b];
                        }
                        for (int b = 0; b < 5; b++) {
                                v[r][b][i] = (tot && tot >= ax.min_bases) ? 100.0 * (double)c[b] / (double)tot : NAN;
                        }
                }
        }
        svg_axis ya = {0.0, 100.0, 0, "% of bases", 0};
        char desc[800], range[96];
        pos_range_text(&ax, range, sizeof(range));
        for (int r = 0; r < nm; r++) {
                struct panel_data pd;
                memset(&pd, 0, sizeof(pd));
                pd.title = mate_title(nm, r);
                pd.x = ax.x;
                pd.y = ya;
                pd.n = ax.nb;
                pd.x0 = ax.x0;
                pd.dx = 1.0;
                pd.nl = 5;
                for (int b = 0; b < 5; b++) {
                        pd.li[b] = (struct lineinfo){base_series[b], report_base_colour[b], report_base_dash[b], 0};
                        pd.v[b] = has_values(v[r][b], ax.nb, 0) ? v[r][b] : NULL;
                }
                snprintf(desc, sizeof(desc),
                         "Percentage of A, C, G, T and N at each position along the read, %s%s. "
                         "One line per base, each with its own dash pattern; the line key is above the figure. "
                         "Positions covered by fewer than " REPORT_MIN_POS_BASES_STR " reads are not drawn.",
                         range, s->is_aligned ? ", aligned part of the read" : "");
                if (figure(o, NULL, "samstat-p4", help_of(4)->title, desc, r, nm, &pd, 1, PH_ONE) != OK) {
                        return FAIL;
                }
        }
        return OK;
}

/* ------------------------------------------------------------------ */
/* plot 5: where along the read do errors occur?                        */
/* ------------------------------------------------------------------ */

static const char* const err_row_title[3] = {"Mismatches", "Insertions", "Deletions"};

#define ERR_MIN_TOP 0.5   /* smallest top of a plot 5 rate axis, in percent */

static double err_rate(const struct collect_pos_bin* b, int row)
{
        uint64_t num = row == 0 ? b->mismatches : row == 1 ? b->insertions : b->deletions;
        uint64_t den = row == 0 ? b->aligned_md : b->aligned;
        return den ? 100.0 * (double)num / (double)den : NAN;
}

static int plot_errors(tld_strbuf* o, const struct stat_collection* s, struct pool* pl)
{
        const int nm = mates(s);
        struct lineinfo li[MAXL];
        struct posaxis ax;
        uint64_t md_total = 0;
        char buf1[48], buf2[48], text[600];

        if (!s->is_aligned) {
                return OK;      /* plots 5 and 6 do not apply to unaligned input */
        }
        const uint64_t mapped = n_mapped(s);
        if (mapped == 0) {
                return na_line(o, help_of(5)->na_no_mapped);
        }
        if (s->total_bases == 0) {
                return na_line(o, help_of(5)->na_no_sequence);
        }
        if (!pos_axis(s, &ax)) {
                return too_few_line(o, s);
        }
        const int nl = get_lines(s, 0, li);   /* mapped classes only */
        for (int r = 0; r < COLLECT_N_MATE; r++) {
                for (int m = 0; m < COLLECT_MAPQ_UNMAPPED; m++) {
                        const struct collect_pos_bin* pb = pos_bins(s, &ax, r, m);
                        for (int i = 0; i < ax.nb; i++) {
                                md_total += pb[i].aligned_md;
                        }
                }
        }
        const int first = md_total ? 0 : 1;          /* the mismatch row is dropped without MD or X */
        const int nrows = 3 - first;

        if (!md_total) {
                if (html_para(o, "plot-na", NULL, report_mismatch_row_na) != OK) {
                        return FAIL;
                }
        } else if (s->reads_without_md > 0) {
                fmt_count(buf1, sizeof(buf1), s->reads_without_md);
                fmt_count(buf2, sizeof(buf2), mapped);
                snprintf(text, sizeof(text), report_note_mismatch_excludes, buf1, buf2);
                if (html_para(o, "plot-note", "Mismatches:", text) != OK) {
                        return FAIL;
                }
        }

        double* v[2][3][MAXL];
        double ymax[3] = {0.0, 0.0, 0.0};
        double pmax[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};   /* per mate and row, for "none detected" */
        int seen[2][3] = {{0, 0, 0}, {0, 0, 0}};                  /* a rate was measured at some position */
        for (int r = 0; r < nm; r++) {
                for (int row = first; row < 3; row++) {
                        for (int k = 0; k < nl; k++) {
                                v[r][row][k] = pool_alloc(pl, (size_t)ax.nb);
                                if (!v[r][row][k]) {
                                        return FAIL;
                                }
                                for (int i = 0; i < ax.nb; i++) {
                                        struct collect_pos_bin sum;
                                        memset(&sum, 0, sizeof(sum));
                                        for (int m = 0; m < COLLECT_MAPQ_UNMAPPED; m++) {
                                                if (li[k].mask & (1u << m)) {
                                                        const struct collect_pos_bin* pb = &pos_bins(s, &ax, r, m)[i];
                                                        sum.mismatches += pb->mismatches;
                                                        sum.insertions += pb->insertions;
                                                        sum.deletions += pb->deletions;
                                                        sum.aligned += pb->aligned;
                                                        sum.aligned_md += pb->aligned_md;
                                                }
                                        }
                                        double x = line_support(s, &ax, r, li[k].mask, i) >= ax.min_bases
                                                           ? err_rate(&sum, row) : NAN;
                                        v[r][row][k][i] = x;
                                        if (isfinite(x)) {
                                                seen[r][row] = 1;
                                                if (x > pmax[r][row]) {
                                                        pmax[r][row] = x;
                                                }
                                                if (x > ymax[row]) {
                                                        ymax[row] = x;
                                                }
                                        }
                                }
                        }
                }
        }
        char names[200], desc[900], range[96];
        char titles[2][3][64];
        pos_range_text(&ax, range, sizeof(range));
        for (int r = 0; r < nm; r++) {
                struct panel_data rows[3];
                int n_none = 0;
                memset(rows, 0, sizeof(rows));
                for (int row = first; row < 3; row++) {
                        struct panel_data* pd = &rows[row - first];
                        snprintf(titles[r][row], sizeof(titles[r][row]), "%s%s (%% of aligned bases)",
                                 nm == 2 ? (r ? "Read 2: " : "Read 1: ") : "", err_row_title[row]);
                        pd->title = titles[r][row];
                        pd->x = ax.x;
                        /* all zero: a whole-number scale (0, 2, 4 %) instead of 0.0-1.0, and a label that says so.
                         * Otherwise the scale reaches at least ERR_MIN_TOP %, so that a handful of events among
                         * thousands of reads (0.004 %) stays a flat line near zero instead of filling the panel. */
                        pd->y = (svg_axis){0.0, ymax[row] > 0.0 ? fmax(svg_nice_ceil(ymax[row] * 1.04, 3), ERR_MIN_TOP) : 4.0, 0, NULL, 0};
                        if (seen[r][row] && pmax[r][row] <= 0.0) {
                                pd->note = "none detected";
                                n_none++;
                        }
                        pd->n = ax.nb;
                        pd->x0 = ax.x0;
                        pd->dx = 1.0;
                        pd->nl = nl;
                        for (int k = 0; k < nl; k++) {
                                pd->li[k] = li[k];
                                pd->v[k] = has_values(v[r][row][k], ax.nb, 0) ? v[r][row][k] : NULL;
                        }
                }
                list_lines(&rows[nrows - 1], names, sizeof(names));
                snprintf(desc, sizeof(desc),
                         "%s per aligned base at each position along the read, %s, aligned part of the read, "
                         "as percent of aligned bases. One panel per error type, sharing the position axis; each panel has its own vertical scale. Lines: %s.%s "
                         "A line is not drawn where fewer than " REPORT_MIN_POS_BASES_STR " of its reads cover the position.",
                         md_total ? "Mismatches, insertions and deletions" : "Insertions and deletions", range, names,
                         n_none ? " A panel labelled none detected has an error rate of zero at every position." : "");
                if (figure(o, s, "samstat-p5", help_of(5)->title, desc, r, nm, rows, nrows, PH_ROW) != OK) {
                        return FAIL;
                }
        }
        return OK;
}

/* ------------------------------------------------------------------ */
/* plot 6: how accurate are the aligned reads?                          */
/* ------------------------------------------------------------------ */

static int plot_identity(tld_strbuf* o, const struct stat_collection* s, struct pool* pl, int shape)
{
        const int nm = mates(s);
        struct lineinfo li[MAXL];
        double* h[2][MAXL];
        const int is_long = (s->layout == COLLECT_LAYOUT_LONG);
        uint64_t all[COLLECT_N_IDENT];

        if (!s->is_aligned) {
                return OK;
        }
        if (shape && (n_mapped(s) == 0 || s->identity_reads == 0)) {
                return OK;
        }
        if (n_mapped(s) == 0) {
                return na_line(o, help_of(6)->na_no_mapped);
        }
        if (s->identity_reads == 0) {
                return na_line(o, help_of(6)->na_no_data);
        }
        const int nl = get_lines(s, 0, li);

        /* 1st percentile of the pooled distribution -> left end of the axis, whole percent, at most 95
         * (the axis is at least 5 points wide: whole-number ticks, and a pile at 100 is a bar, not the frame) */
        memset(all, 0, sizeof(all));
        uint64_t total = 0;
        for (int r = 0; r < COLLECT_N_MATE; r++) {
                for (int m = 0; m < COLLECT_MAPQ_UNMAPPED; m++) {
                        for (int b = 0; b < COLLECT_N_IDENT; b++) {
                                all[b] += s->acc->identity[r][m][b];
                                total += s->acc->identity[r][m][b];
                        }
                }
        }
        if (total == 0) {
                return na_line(o, help_of(6)->na_no_data);
        }
        int p1 = 0;
        {
                uint64_t cum = 0;
                for (int b = 0; b < COLLECT_N_IDENT; b++) {
                        cum += all[b];
                        if (cum * 100 >= total) {
                                p1 = b;
                                break;
                        }
                }
        }
        int lo_pct = p1 / 10;
        if (lo_pct > 95) {
                lo_pct = 95;
        }

        /*
         * Short layout: 1 % bins (a 150 bp read moves in steps of 0.67 %, finer bins would alias).
         * Long layout: the narrowest of 0.1, 0.2, 0.5 and 1 % that gives at most 50 bins over the
         * axis, so that a few thousand reads still draw a curve, not a comb. Each bin is a sum of
         * whole stored 0.1 % bins (k of them). The axis ends at 100 and the stored exact-100 % bin
         * joins the last bin, [100 - bin width, 100].
         */
        int k10 = 10;
        if (is_long) {
                static const int widths[4] = {1, 2, 5, 10};
                for (int w = 0; w < 4; w++) {
                        k10 = widths[w];
                        if ((1000 - lo_pct * 10) / k10 <= 50) {
                                break;
                        }
                }
        }
        const double bw = k10 / 10.0;
        const int n = (1000 - lo_pct * 10) / k10;
        double ymax = 0.0;
        for (int r = 0; r < nm; r++) {
                for (int k = 0; k < nl; k++) {
                        h[r][k] = pool_alloc(pl, (size_t)n);
                        if (!h[r][k]) {
                                return FAIL;
                        }
                        for (int m = 0; m < COLLECT_MAPQ_UNMAPPED; m++) {
                                if (!(li[k].mask & (1u << m))) {
                                        continue;
                                }
                                for (int b = lo_pct * 10; b < COLLECT_N_IDENT; b++) {
                                        uint64_t c = s->acc->identity[r][m][b];
                                        if (!c) {
                                                continue;
                                        }
                                        int d = (b - lo_pct * 10) / k10;
                                        if (d > n - 1) {
                                                d = n - 1;      /* exactly 100 %: the last bin */
                                        }
                                        h[r][k][d] += (double)c;
                                }
                        }
                        for (int d = 0; d < n; d++) {
                                ymax = h[r][k][d] > ymax ? h[r][k][d] : ymax;
                        }
                }
        }
        if (shape && !shape_scale(h, nm, nl, 0, n, &ymax)) {
                return OK;
        }
        svg_axis xa = {(double)lo_pct, 100.0, 0, "identity to the reference (%)", 1};
        svg_axis ya = shape ? share_axis(ymax) : (svg_axis){0.0, count_top(ymax), 0, "number of reads", 1};
        char names[200], desc[900], med[48];

        for (int r = 0; r < nm; r++) {
                struct panel_data pd;
                memset(&pd, 0, sizeof(pd));
                pd.title = mate_title(nm, r);
                pd.x = xa;
                pd.y = ya;
                pd.step = 1;
                pd.n = n;
                pd.x0 = (double)lo_pct;
                pd.dx = bw;
                pd.nl = nl;
                for (int k = 0; k < nl; k++) {
                        pd.li[k] = li[k];
                        pd.v[k] = h[r][k] && has_values(h[r][k], n, 1) ? h[r][k] : NULL;
                }
                if (!shape && isfinite(s->median_identity)) {
                        fmt_dec(med, sizeof(med), 1, s->median_identity);
                        pd.has_vline = 1;
                        pd.vline = s->median_identity;
                        snprintf(pd.vlabel, sizeof(pd.vlabel), "median %s%%", med);
                } else {
                        snprintf(med, sizeof(med), "unavailable");
                }
                list_lines(&pd, names, sizeof(names));
                snprintf(desc, sizeof(desc),
                         "Number of aligned reads by identity to the reference, from %d%% to 100%%, in bins of %s. "
                         "Lines: %s. The dashed vertical line marks the median identity, %s%%. "
                         "The last bin includes the reads with exactly 100%% identity.",
                         lo_pct, k10 == 10 ? "1 percentage point" : k10 == 5 ? "0.5 percentage points"
                                 : k10 == 2 ? "0.2 percentage points" : "0.1 percentage points", names, med);
                if (shape) {
                        shape_desc(desc, sizeof(desc), "by identity to the reference", names);
                }
                if (figure(o, s, shape ? "samstat-p6s" : "samstat-p6", help_of(6)->title, desc, r, nm, &pd, 1, PH_ONE) != OK) {
                        return FAIL;
                }
        }
        return OK;
}

/* ------------------------------------------------------------------ */
/* public interface                                                     */
/* ------------------------------------------------------------------ */

static int render(tld_strbuf* o, const struct stat_collection* s, int plot_id, int shape)
{
        struct pool pl = {NULL, 0, 0};
        int rc = FAIL;

        if (!o || !s || plot_id < 1 || plot_id > REPORT_N_PLOTS) {
                ERROR_MSG("report_plot_render: bad arguments (plot %d)", plot_id);
        }
        if (s->n_total == 0) {
                /* an empty input is valid; plots 5 and 6 are absent on unaligned input */
                if (plot_id >= 5 && !s->is_aligned) {
                        return OK;
                }
                return na_line(o, NA_NO_READS);
        }
        pl.cap = (size_t)(2 * 3 * MAXL + 4) * MAXB;
        pl.d = calloc(pl.cap, sizeof(double));
        if (!pl.d) {
                ERROR_MSG("report_plot_render: out of memory");
        }
        switch (plot_id) {
        case 1:
                rc = plot_length(o, s, &pl, shape);
                break;
        case 2:
                rc = plot_readq(o, s, &pl, shape);
                break;
        case 3:
                rc = plot_posq(o, s, &pl);
                break;
        case 4:
                rc = plot_comp(o, s, &pl);
                break;
        case 5:
                rc = plot_errors(o, s, &pl);
                break;
        default:
                rc = plot_identity(o, s, &pl, shape);
                break;
        }
        free(pl.d);
        return rc;
ERROR:
        free(pl.d);
        return FAIL;
}

int report_plot_render(tld_strbuf* o, const struct stat_collection* s, int plot_id)
{
        return render(o, s, plot_id, 0);
}

int report_plot_render_shape(tld_strbuf* o, const struct stat_collection* s, int plot_id)
{
        if (!o || !s || (plot_id != 1 && plot_id != 2 && plot_id != 6)) {
                return OK;
        }
        if (s->n_total == 0 || !s->is_aligned) {
                return OK;
        }
        return render(o, s, plot_id, 1);
}

int report_mapq_swatch(tld_strbuf* o, int mapq_class)
{
        static const char* const colour_name[4] = {"blue", "orange", "vermillion", "grey"};
        static const char* const dash_name[4] = {"solid", "dashed", "dotted", "dash-dot"};
        char id[32], title[64], desc[96];

        if (!o || mapq_class < 0 || mapq_class >= COLLECT_N_MAPQ) {
                ERROR_MSG("report_mapq_swatch: bad arguments");
        }
        /* The class name is the row header next to the swatch: the figure names only the style */
        snprintf(id, sizeof(id), "samstat-mapq-%d", mapq_class);
        snprintf(title, sizeof(title), "%s %s line", dash_name[mapq_class], colour_name[mapq_class]);
        snprintf(desc, sizeof(desc), "Line sample: %s, %s.", colour_name[mapq_class], dash_name[mapq_class]);
        return svg_swatch(o, id, report_mapq_colour[mapq_class], report_mapq_dash[mapq_class], title, desc);
ERROR:
        return FAIL;
}
