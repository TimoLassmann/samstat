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

/* One panel is FIG_W units wide. A figure is shown at its natural size when the container is at
 * least that wide and shrinks with its text below that: at a 360 px screen (328 px of content)
 * the 13 unit text is about 11.8 px. */
#define FIG_W 420.0
#define FONT 14.0
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

/* Lengths as in the summary strip: exact below 1 kb, "~" above (the value is a bin edge, within
 * 1%), bp below 10 kb, kb above. */
static void fmt_len(char* buf, size_t n, uint64_t v)
{
        char t[48];
        if (v < 1000) {
                fmt_count(t, sizeof(t), v);
                snprintf(buf, n, "%s bp", t);
        } else if (v < 10000) {
                fmt_count(t, sizeof(t), v);
                snprintf(buf, n, "~%s bp", t);
        } else {
                char d[32];
                fmt_dec(d, sizeof(d), 1, (double)v / 1000.0);
                /* thousands separators in the integer part of the kb value */
                char* dot = strchr(d, '.');
                size_t il = dot ? (size_t)(dot - d) : strlen(d);
                char ip[32];
                snprintf(ip, sizeof(ip), "%.*s", (int)il, d);
                fmt_count(t, sizeof(t), (uint64_t)strtoull(ip, NULL, 10));
                snprintf(buf, n, "~%s%s kb", t, dot ? dot : "");
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
};

/* Top of a count axis: a little headroom over the tallest bar, and never below 4 so that a handful of
 * reads does not get ticks at 0.2 reads. */
static double count_top(double ymax)
{
        double t = svg_nice_ceil(ymax * 1.04, 5);
        return t < 4.0 ? 4.0 : t;
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

static int figure(tld_strbuf* o, const char* id_base, const char* title, const char* desc,
                  int mate, int nm, const struct panel_data* rows, int nrows, double ph)
{
        svg_grid g;
        svg_panel p;
        char id[64];
        char ttl[256];
        char sertitle[128];
        const int titles = 1;

        svg_grid_init(&g, nrows, 1, FIG_W, 300.0, FONT, titles);
        double h = g.mt + g.mb + (nrows - 1) * g.gap_y + nrows * ph;
        svg_grid_init(&g, nrows, 1, FIG_W, h, FONT, titles);

        if (nm == 2) {
                snprintf(id, sizeof(id), "%s-r%d", id_base, mate + 1);
                snprintf(ttl, sizeof(ttl), "%s Read %d", title, mate + 1);
        } else {
                snprintf(id, sizeof(id), "%s", id_base);
                snprintf(ttl, sizeof(ttl), "%s", title);
        }
        TRY(svg_begin(o, id, FIG_W, h, FONT, ttl, desc));
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
        svg_axis x;
};

/*
 * x axis of plots 3-5 (3.2). Short layout: absolute, from position 1 to the last bin with any
 * base (over all classes and both mates). Long layout: relative, 0-100 %. Returns 0 when no
 * base was seen.
 */
static int pos_axis(const struct stat_collection* s, struct posaxis* a)
{
        const struct collect_acc* acc = s->acc;
        static char lab_abs[64], lab_rel[64];

        a->rel = (s->layout == COLLECT_LAYOUT_LONG);
        snprintf(lab_abs, sizeof(lab_abs), "position in read (bp)%s", s->is_aligned ? ", aligned part" : "");
        snprintf(lab_rel, sizeof(lab_rel), "position in read (%% of length)%s", s->is_aligned ? ", aligned part" : "");
        if (a->rel) {
                int any = 0;
                for (int r = 0; r < COLLECT_N_MATE && !any; r++) {
                        for (int m = 0; m < COLLECT_N_MAPQ && !any; m++) {
                                for (int i = 0; i < COLLECT_REL_BINS && !any; i++) {
                                        const struct collect_pos_bin* b = &acc->pos_rel[r][m][i];
                                        for (int c = 0; c < COLLECT_N_BASE; c++) {
                                                any |= b->comp[c] != 0;
                                        }
                                }
                        }
                }
                a->nb = COLLECT_REL_BINS;
                a->x0 = 0.5;
                a->x = (svg_axis){0.0, 100.0, 0, lab_rel};
                return any;
        }
        int last = -1;
        for (int r = 0; r < COLLECT_N_MATE; r++) {
                for (int m = 0; m < COLLECT_N_MAPQ; m++) {
                        for (int i = COLLECT_ABS_BINS - 1; i > last; i--) {
                                const struct collect_pos_bin* b = &acc->pos_abs[r][m][i];
                                uint64_t t = 0;
                                for (int c = 0; c < COLLECT_N_BASE; c++) {
                                        t += b->comp[c];
                                }
                                if (t) {
                                        last = i;
                                        break;
                                }
                        }
                }
        }
        a->nb = last + 1;
        a->x0 = 1.0;
        /* positions 1..nb are drawn at x = 1..nb; half a bin of room on each side keeps the lines off the
         * frame and lets the ticks fall on 50, 100, 150 */
        a->x = (svg_axis){0.5, (double)(a->nb > 0 ? a->nb : 1) + 0.5, 0, lab_abs};
        return a->nb > 0;
}

static const struct collect_pos_bin* pos_bins(const struct stat_collection* s, const struct posaxis* a, int r, int m)
{
        return a->rel ? s->acc->pos_rel[r][m] : s->acc->pos_abs[r][m];
}

/* ------------------------------------------------------------------ */
/* plot 1: how long are the reads?                                      */
/* ------------------------------------------------------------------ */

static double len_edge(int is_long, int b)
{
        return (double)(is_long ? collect_display_edge(b) : collect_length_edge(b));
}

/* Upper edge of bin b: the next lower edge; the overflow bin gets one more grid step (1%). */
static double len_upper(int is_long, int b)
{
        int nb = is_long ? collect_n_display_bins() : COLLECT_N_LEN_BINS;
        if (b + 1 < nb) {
                return len_edge(is_long, b + 1);
        }
        return len_edge(is_long, b) * 1.01 + 1.0;
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

static int plot_length(tld_strbuf* o, const struct stat_collection* s, struct pool* pl)
{
        const int is_long = (s->layout == COLLECT_LAYOUT_LONG);
        const int nb = is_long ? collect_n_display_bins() : COLLECT_N_LEN_BINS;
        const int nm = mates(s);
        struct lineinfo li[MAXL];
        const int nl = get_lines(s, 1, li);
        double* cnt[2][MAXL];
        int bmin = nb, bmax = -1;
        double ymax = 0.0;

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
                                                cnt[r][k][is_long ? collect_display_bin(b) : b] += (double)v;
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
                return na_line(o, NA_NO_READS);
        }

        const int n = bmax - bmin + 1;
        double* xe = pool_alloc(pl, (size_t)n + 1);
        if (!xe) {
                return FAIL;
        }
        for (int i = 0; i < n; i++) {
                xe[i] = len_edge(is_long, bmin + i);
        }
        xe[n] = len_upper(is_long, bmax);

        svg_axis xa;
        if (is_long) {
                double lo = fmax(xe[0], 1.0);
                double hi = xe[n];
                if (hi < lo * 10.0) {   /* at least one decade, centred on the data */
                        double mid = sqrt(lo * hi);
                        lo = mid / sqrt(10.0);
                        hi = mid * sqrt(10.0);
                }
                xa = (svg_axis){nice125_down(fmax(lo, 1.0)), nice125_up(hi), 1, "read length (bp), log scale"};
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
                xa = (svg_axis){lo, hi, 0, "read length (bp)"};
        }
        svg_axis ya = {0.0, count_top(ymax), 0, "number of reads"};

        char lo_s[48], hi_s[48], n50_s[48], names[200], desc[900];
        fmt_len(lo_s, sizeof(lo_s), s->min_length);
        fmt_len(hi_s, sizeof(hi_s), s->max_length);
        fmt_len(n50_s, sizeof(n50_s), s->n50);

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
                        pd.v[k] = has_values(cnt[r][k] + bmin, n, 1) ? cnt[r][k] + bmin : NULL;
                }
                if (s->n50 > 0) {
                        pd.has_vline = 1;
                        pd.vline = (double)s->n50;
                        snprintf(pd.vlabel, sizeof(pd.vlabel), "N50 %s", n50_s);
                }
                list_lines(&pd, names, sizeof(names));
                snprintf(desc, sizeof(desc),
                         "Number of reads per read-length bin%s. Read lengths range from %s to %s. Lines: %s. "
                         "The dashed vertical line marks N50, %s: half of all bases are in reads at least this long.",
                         is_long ? " on a log scale" : "", lo_s, hi_s, names, n50_s);
                if (figure(o, "samstat-p1", help_of(1)->title, desc, r, nm, &pd, 1, PH_ONE) != OK) {
                        return FAIL;
                }
        }
        return OK;
}

/* ------------------------------------------------------------------ */
/* plot 2: how good is each read?                                       */
/* ------------------------------------------------------------------ */

static int plot_readq(tld_strbuf* o, const struct stat_collection* s, struct pool* pl)
{
        const int nm = mates(s);
        struct lineinfo li[MAXL];
        const int nl = get_lines(s, 1, li);
        double* cnt[2][MAXL];
        int qmax = -1;
        double ymax = 0.0;

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
                                        ymax = cnt[r][k][q] > ymax ? cnt[r][k][q] : ymax;
                                }
                        }
                }
        }
        if (qmax < 0) {
                return na_line(o, help_of(2)->na_no_quality);
        }
        const int n = qmax + 1;
        svg_axis xa = {0.0, (double)n, 0, "mean base quality of the read (Phred Q)"};
        svg_axis ya = {0.0, count_top(ymax), 0, "number of reads"};
        char names[200], desc[700];

        for (int r = 0; r < nm; r++) {
                struct panel_data pd;
                memset(&pd, 0, sizeof(pd));
                pd.title = mate_title(nm, r);
                pd.x = xa;
                pd.y = ya;
                pd.step = 1;
                pd.n = n;
                pd.x0 = 0.0;
                pd.dx = 1.0;
                pd.nl = nl;
                for (int k = 0; k < nl; k++) {
                        pd.li[k] = li[k];
                        pd.v[k] = has_values(cnt[r][k], n, 1) ? cnt[r][k] : NULL;
                }
                list_lines(&pd, names, sizeof(names));
                snprintf(desc, sizeof(desc),
                         "Number of reads by their mean base quality, from Q0 to Q%d. Lines: %s. "
                         "Qualities are averaged as error probabilities.", n, names);
                if (figure(o, "samstat-p2", help_of(2)->title, desc, r, nm, &pd, 1, PH_ONE) != OK) {
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
                return na_line(o, NA_NO_READS);
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
                                double q = collect_mean_q(err, bases);   /* NaN when there is no base */
                                v[r][k][i] = q;
                                if (isfinite(q) && q > qmax) {
                                        qmax = q;
                                }
                        }
                }
        }
        svg_axis ya = {0.0, fmax(10.0, ceil(qmax / 10.0 - 1e-9) * 10.0), 0, "mean base quality (Phred Q)"};
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
                         "Gaps mean that no base with a quality was seen.",
                         range, s->is_aligned ? ", aligned part of the read" : "", names);
                if (figure(o, "samstat-p3", help_of(3)->title, desc, r, nm, &pd, 1, PH_ONE) != OK) {
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
static const char* const base_colour_name[5] = {"green", "reddish purple", "black", "brown", "sky blue"};

static int plot_comp(tld_strbuf* o, const struct stat_collection* s, struct pool* pl)
{
        const int nm = mates(s);
        struct posaxis ax;
        double* v[2][5];

        if (!pos_axis(s, &ax)) {
                return na_line(o, NA_NO_READS);
        }
        /* the key: one row of five swatches above the panels */
        {
                svg_swatch_item items[5];
                char desc[300];
                int l = snprintf(desc, sizeof(desc), "Line key for the bases:");
                for (int b = 0; b < 5; b++) {
                        items[b] = (svg_swatch_item){base_name[b], report_base_colour[b], report_base_dash[b]};
                        l += snprintf(desc + l, sizeof(desc) - (size_t)l, "%s %s %s, %s", b ? ";" : "", base_name[b],
                                      base_colour_name[b], report_base_dash[b][0] ? "dashed line" : "solid line");
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
                                v[r][b][i] = tot ? 100.0 * (double)c[b] / (double)tot : NAN;
                        }
                }
        }
        svg_axis ya = {0.0, 100.0, 0, "% of bases"};
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
                         "One line per base; the line key is above the figure.",
                         range, s->is_aligned ? ", aligned part of the read" : "");
                if (figure(o, "samstat-p4", help_of(4)->title, desc, r, nm, &pd, 1, PH_ONE) != OK) {
                        return FAIL;
                }
        }
        return OK;
}

/* ------------------------------------------------------------------ */
/* plot 5: where along the read do errors occur?                        */
/* ------------------------------------------------------------------ */

static const char* const err_row_title[3] = {"Mismatches", "Insertions", "Deletions"};

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
        if (!pos_axis(s, &ax)) {
                return na_line(o, NA_NO_READS);
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
                                        double x = err_rate(&sum, row);
                                        v[r][row][k][i] = x;
                                        if (isfinite(x) && x > ymax[row]) {
                                                ymax[row] = x;
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
                memset(rows, 0, sizeof(rows));
                for (int row = first; row < 3; row++) {
                        struct panel_data* pd = &rows[row - first];
                        snprintf(titles[r][row], sizeof(titles[r][row]), "%s%s (%% of aligned bases)",
                                 nm == 2 ? (r ? "Read 2: " : "Read 1: ") : "", err_row_title[row]);
                        pd->title = titles[r][row];
                        pd->x = ax.x;
                        pd->y = (svg_axis){0.0, svg_nice_ceil(ymax[row] * 1.04, 3), 0, NULL};
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
                         "as percent of aligned bases. One panel per error type, sharing the position axis. Lines: %s.",
                         md_total ? "Mismatches, insertions and deletions" : "Insertions and deletions", range, names);
                if (figure(o, "samstat-p5", help_of(5)->title, desc, r, nm, rows, nrows, PH_ROW) != OK) {
                        return FAIL;
                }
        }
        return OK;
}

/* ------------------------------------------------------------------ */
/* plot 6: how accurate are the aligned reads?                          */
/* ------------------------------------------------------------------ */

static int plot_identity(tld_strbuf* o, const struct stat_collection* s, struct pool* pl)
{
        const int nm = mates(s);
        struct lineinfo li[MAXL];
        double* h[2][MAXL];
        const int is_long = (s->layout == COLLECT_LAYOUT_LONG);
        uint64_t all[COLLECT_N_IDENT];

        if (!s->is_aligned) {
                return OK;
        }
        if (n_mapped(s) == 0) {
                return na_line(o, help_of(6)->na_no_mapped);
        }
        if (s->identity_reads == 0) {
                return na_line(o, help_of(6)->na_no_data);
        }
        const int nl = get_lines(s, 0, li);

        /* 1st percentile of the pooled distribution -> left end of the axis, whole percent, at most 99 */
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
        if (lo_pct > 99) {
                lo_pct = 99;
        }

        /*
         * Long layout: the stored 0.1 % bins. Short layout: 1 % bins (ten stored bins summed); the
         * exact-100 % bin stays a bin of its own. Either way the last bin is drawn just right of 100.
         */
        const double bw = is_long ? 0.1 : 1.0;
        const int n = is_long ? COLLECT_N_IDENT - lo_pct * 10 : (100 - lo_pct) + 1;
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
                                        int d = is_long ? b - lo_pct * 10 : (b == 1000 ? n - 1 : b / 10 - lo_pct);
                                        h[r][k][d] += (double)c;
                                }
                        }
                        for (int d = 0; d < n; d++) {
                                ymax = h[r][k][d] > ymax ? h[r][k][d] : ymax;
                        }
                }
        }
        svg_axis xa = {(double)lo_pct, 100.0 + bw, 0, "identity to the reference (%)"};
        svg_axis ya = {0.0, count_top(ymax), 0, "number of reads"};
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
                pd.x_tick_max = 100.0;
                pd.nl = nl;
                for (int k = 0; k < nl; k++) {
                        pd.li[k] = li[k];
                        pd.v[k] = has_values(h[r][k], n, 1) ? h[r][k] : NULL;
                }
                if (isfinite(s->median_identity)) {
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
                         "The last bin holds the reads with exactly 100%% identity.",
                         lo_pct, is_long ? "0.1 percentage points" : "1 percentage point", names, med);
                if (figure(o, "samstat-p6", help_of(6)->title, desc, r, nm, &pd, 1, PH_ONE) != OK) {
                        return FAIL;
                }
        }
        return OK;
}

/* ------------------------------------------------------------------ */
/* public interface                                                     */
/* ------------------------------------------------------------------ */

int report_plot_render(tld_strbuf* o, const struct stat_collection* s, int plot_id)
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
                rc = plot_length(o, s, &pl);
                break;
        case 2:
                rc = plot_readq(o, s, &pl);
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
                rc = plot_identity(o, s, &pl);
                break;
        }
        free(pl.d);
        return rc;
ERROR:
        free(pl.d);
        return FAIL;
}

int report_mapq_swatch(tld_strbuf* o, int mapq_class)
{
        static const char* const colour_name[4] = {"blue", "orange", "vermillion", "grey"};
        static const char* const dash_name[4] = {"solid", "dashed", "dotted", "dash-dot"};
        char id[32], desc[96];

        if (!o || mapq_class < 0 || mapq_class >= COLLECT_N_MAPQ) {
                ERROR_MSG("report_mapq_swatch: bad arguments");
        }
        snprintf(id, sizeof(id), "samstat-mapq-%d", mapq_class);
        snprintf(desc, sizeof(desc), "Line sample: %s, %s.", colour_name[mapq_class], dash_name[mapq_class]);
        return svg_swatch(o, id, report_mapq_colour[mapq_class], report_mapq_dash[mapq_class],
                          report_mapq_help[mapq_class].label, desc);
ERROR:
        return FAIL;
}
