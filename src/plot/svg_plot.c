#include "svg_plot.h"

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

/* Neutral chrome. Series colours come from the caller. */
#define SVG_TEXT_COLOR "#1a1a1a"   /* 16:1 on white */
#define SVG_AXIS_COLOR "#595959"   /*  7:1 on white */
#define SVG_GRID_COLOR "#e4e4e4"   /* decorative */
#define SVG_FONT_STACK "system-ui, -apple-system, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif"
#define SVG_LINE_WIDTH 2.0
#define SVG_MIN_BAR 6.0           /* an isolated histogram bin is drawn at least this wide (user units) */

/*
 * Like tld's RUN, but silent: RUN prints the stringified expression through a
 * printf-style call, and these expressions contain '%' (format strings).
 */
#define TRY(EXP) do { if ((EXP) != OK) { goto ERROR; } } while (0)

/* ------------------------------------------------------------------ */
/* output helpers                                                       */
/* ------------------------------------------------------------------ */

/* Locale independent, at most one decimal, no trailing zeros, never "-0" or "nan". */
static const char* num(char* buf, double v)
{
        if (!isfinite(v)) {
                v = 0.0;
        }
        if (v > 1e6) {
                v = 1e6;
        } else if (v < -1e6) {
                v = -1e6;
        }
        snprintf(buf, 32, "%.1f", v);
        char* p = buf;
        for (; *p; p++) {
                if (*p == ',') {
                        *p = '.';
                }
        }
        if (strchr(buf, '.')) {
                char* e = buf + strlen(buf) - 1;
                while (*e == '0') {
                        *e-- = 0;
                }
                if (*e == '.') {
                        *e = 0;
                }
        }
        if (strcmp(buf, "-0") == 0) {
                strcpy(buf, "0");
        }
        return buf;
}

static int put_n(tld_strbuf* o, const char* s, size_t n)
{
        char buf[256];
        while (n) {
                size_t k = n > sizeof(buf) - 1 ? sizeof(buf) - 1 : n;
                memcpy(buf, s, k);
                buf[k] = 0;
                TRY(tld_append(o, buf));
                s += k;
                n -= k;
        }
        return OK;
ERROR:
        return FAIL;
}

/* XML text/attribute escaping; control characters become spaces. */
static int put_escaped(tld_strbuf* o, const char* s)
{
        const char* run = s;
        if (!s) {
                return OK;
        }
        for (const char* c = s;; c++) {
                const char* rep = NULL;
                switch (*c) {
                case 0:
                        return put_n(o, run, (size_t)(c - run));
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
                case '\'':
                        rep = "&#39;";
                        break;
                default:
                        if ((unsigned char)*c < 0x20) {
                                rep = " ";
                        }
                        break;
                }
                if (rep) {
                        TRY(put_n(o, run, (size_t)(c - run)));
                        TRY(put_n(o, rep, strlen(rep)));
                        run = c + 1;
                }
        }
ERROR:
        return FAIL;
}

/*
 * Tiny formatter: %f = double as a coordinate, %s = raw string, %e = escaped
 * string, %d = int, %% = percent sign. Everything else is copied verbatim.
 */
static int pf(tld_strbuf* o, const char* fmt, ...)
{
        va_list ap;
        char tmp[40];
        const char* run = fmt;

        va_start(ap, fmt);
        for (const char* c = fmt;; c++) {
                if (*c && *c != '%') {
                        continue;
                }
                TRY(put_n(o, run, (size_t)(c - run)));
                if (!*c) {
                        break;
                }
                c++;
                switch (*c) {
                case 'f':
                        TRY(tld_append(o, (char*)num(tmp, va_arg(ap, double))));
                        break;
                case 's': {
                        const char* s = va_arg(ap, const char*);
                        TRY(put_n(o, s ? s : "", s ? strlen(s) : 0));
                        break;
                }
                case 'e':
                        TRY(put_escaped(o, va_arg(ap, const char*)));
                        break;
                case 'd':
                        snprintf(tmp, sizeof(tmp), "%d", va_arg(ap, int));
                        TRY(tld_append(o, tmp));
                        break;
                case '%':
                        TRY(tld_append_char(o, '%'));
                        break;
                default:
                        ERROR_MSG("svg_plot: bad format");
                }
                run = c + 1;
        }
        va_end(ap);
        return OK;
ERROR:
        va_end(ap);
        return FAIL;
}

/* ------------------------------------------------------------------ */
/* ticks                                                                */
/* ------------------------------------------------------------------ */

double svg_nice_step(double range, int target)
{
        if (!(range > 0.0) || !isfinite(range)) {
                return 1.0;
        }
        if (target < 1) {
                target = 1;
        }
        double raw = range / target;
        double mag = pow(10.0, floor(log10(raw)));
        double n = raw / mag;
        double m = n < 1.5 ? 1.0 : n < 3.0 ? 2.0 : n < 7.0 ? 5.0 : 10.0;
        return m * mag;
}

double svg_nice_ceil(double max, int target)
{
        if (!(max > 0.0) || !isfinite(max)) {
                return 1.0;
        }
        double step = svg_nice_step(max, target);
        return ceil(max / step - 1e-9) * step;
}

/* A usable copy of an axis: finite, hi > lo, and lo > 0 on a log axis. */
static svg_axis sane(const svg_axis* in)
{
        svg_axis a = *in;
        if (!isfinite(a.lo) || !isfinite(a.hi)) {
                a.lo = 0.0;
                a.hi = 1.0;
        }
        if (a.log) {
                if (!(a.lo > 0.0)) {
                        a.lo = 1.0;
                }
                if (!(a.hi > a.lo)) {
                        a.hi = a.lo * 10.0;
                }
        } else if (!(a.hi > a.lo)) {
                a.hi = a.lo + 1.0;
        }
        return a;
}

int svg_ticks(const svg_axis* in, int target, double* v, int* major, int max, double* step_out)
{
        if (!in || !v || !major || max < 1) {
                return -1;
        }
        svg_axis a = sane(in);
        int n = 0;
        if (target < 2) {
                target = 2;
        } else if (target > 40) {
                target = 40;
        }

        if (!a.log) {
                double step = svg_nice_step(a.hi - a.lo, target);
                double k0 = ceil(a.lo / step - 1e-9);
                for (double k = k0; k * step <= a.hi + step * 1e-9 && n < max; k += 1.0) {
                        v[n] = fabs(k) < 0.5 ? 0.0 : k * step;
                        major[n++] = 1;
                }
                if (step_out) {
                        *step_out = step;
                }
                return n;
        }

        /*
         * Log axis. With room: 1, 2, 5 of every decade labelled. Less room:
         * decades labelled, 2 and 5 as short unlabelled ticks while the axis
         * spans at most 6 decades. Less still: every k-th decade.
         */
        double span = log10(a.hi) - log10(a.lo);
        int mode = span * 3.0 <= target ? 0 : span <= target ? 1 : 2;
        int stride = mode == 2 ? (int)ceil(span / target) : 1;
        int e0 = (int)floor(log10(a.lo) + 1e-9);
        int e1 = (int)ceil(log10(a.hi) - 1e-9);
        static const double mant[3] = {1.0, 2.0, 5.0};

        if (step_out) {
                *step_out = 0.0;
        }
        /* A short axis (a lone bin: one or two decades) with fewer than 3 labelled ticks is hard to
         * read: label the 2 and the 5 as well. */
        for (int pass = 0; pass < 2; pass++) {
                int labelled = 0;
                n = 0;
                for (int e = e0; e <= e1 && n < max; e++) {
                        for (int m = 0; m < 3 && n < max; m++) {
                                double val = mant[m] * pow(10.0, e);
                                if (val < a.lo * (1.0 - 1e-9) || val > a.hi * (1.0 + 1e-9)) {
                                        continue;
                                }
                                if (mode == 0) {
                                        major[n] = 1;
                                } else if (m == 0) {
                                        if (mode == 2 && ((e % stride) + stride) % stride != 0) {
                                                continue;
                                        }
                                        major[n] = 1;
                                } else if (mode == 1 && span <= 6.0) {
                                        major[n] = 0;
                                } else {
                                        continue;
                                }
                                labelled += major[n];
                                v[n++] = val;
                        }
                }
                if (pass == 0 && mode == 1 && span <= 2.0 && labelled < 3) {
                        mode = 0;
                } else {
                        break;
                }
        }
        return n;
}

void svg_tick_label(const svg_axis* in, double step, double v, char* buf, size_t n)
{
        svg_axis a = sane(in);
        const char* suffix = "";
        double unit = 1.0;
        int decimals = 0;

        if (fabs(v) < 1e-12) {
                v = 0.0;
        }
        if (a.log) {
                double av = fabs(v);
                if (av >= 1e9) {
                        unit = 1e9;
                        suffix = " G";
                } else if (av >= 1e6) {
                        unit = 1e6;
                        suffix = " M";
                } else if (av >= 1e3) {
                        unit = 1e3;
                        suffix = " k";
                }
                snprintf(buf, n, "%g%s", v / unit, suffix);
        } else {
                double amax = fmax(fabs(a.lo), fabs(a.hi));
                if (amax >= 1e6) {
                        unit = 1e6;
                        suffix = " M";
                } else if (amax >= 1e4) {
                        unit = 1e3;
                        suffix = " k";
                }
                double s = step / unit;
                if (s > 0.0 && s < 1.0 - 1e-9) {
                        decimals = (int)ceil(-log10(s) - 1e-9);
                        if (decimals > 4) {
                                decimals = 4;
                        }
                }
                double scaled = v / unit;
                if (fabs(scaled) < 0.5 * pow(10.0, -decimals)) {
                        scaled = 0.0;
                }
                snprintf(buf, n, "%.*f%s", decimals, scaled, scaled == 0.0 ? "" : suffix);
        }
        for (char* c = buf; *c; c++) {
                if (*c == ',') {
                        *c = '.';
                }
        }
}

/* ------------------------------------------------------------------ */
/* coordinates                                                          */
/* ------------------------------------------------------------------ */

static double frac(const svg_axis* a, double v)
{
        double f;
        if (!(v == v)) {
                return 0.0;
        }
        if (a->log) {
                f = v > a->lo ? (log10(v) - log10(a->lo)) / (log10(a->hi) - log10(a->lo)) : 0.0;
        } else {
                f = (v - a->lo) / (a->hi - a->lo);
        }
        return f < 0.0 ? 0.0 : f > 1.0 ? 1.0 : f;
}

static double px(const svg_panel* p, double v)
{
        return p->left + frac(&p->x, v) * p->w;
}

static double py(const svg_panel* p, double v)
{
        return p->top + (1.0 - frac(&p->y, v)) * p->h;
}

static double series_x(const svg_series* s, int i)
{
        return s->x ? s->x[i] : s->x0 + i * s->dx;
}

/* ------------------------------------------------------------------ */
/* document                                                             */
/* ------------------------------------------------------------------ */

static int begin_ex(tld_strbuf* o, const char* id, double w, double h, double font,
                    const char* title, const char* desc, const char* style)
{
        if (!id || !*id || !(w > 0.0) || !(h > 0.0) || !(font >= SVG_DEFAULT_FONT)) {
                ERROR_MSG("svg_plot: bad svg parameters");
        }
        TRY(pf(o, "<svg xmlns=\"http://www.w3.org/2000/svg\" role=\"img\" aria-labelledby=\"%e-t %e-d\" "
               "width=\"%f\" height=\"%f\" viewBox=\"0 0 %f %f\" style=\"%s\" "
               "font-family=\"%s\" font-size=\"%f\" fill=\"" SVG_TEXT_COLOR "\">",
               id, id, w, h, w, h, style, SVG_FONT_STACK, font));
        TRY(pf(o, "<title id=\"%e-t\">%e</title><desc id=\"%e-d\">%e</desc>", id, title, id, desc));
        return OK;
ERROR:
        return FAIL;
}

int svg_begin(tld_strbuf* o, const char* id, double w, double h, double font,
              const char* title, const char* desc)
{
        TRY(begin_ex(o, id, w, h, font, title, desc, "max-width:100%;height:auto;vertical-align:top"));
        TRY(pf(o, "<rect width=\"%f\" height=\"%f\" fill=\"#fff\"/>", w, h));
        return OK;
ERROR:
        return FAIL;
}

int svg_end(tld_strbuf* o)
{
        return pf(o, "</svg>\n");
}

/* ------------------------------------------------------------------ */
/* grid and panels                                                      */
/* ------------------------------------------------------------------ */

void svg_grid_init(svg_grid* g, int rows, int cols, double w, double h, double font, int has_titles)
{
        g->rows = rows < 1 ? 1 : rows;
        g->cols = cols < 1 ? 1 : cols;
        g->width = w;
        g->height = h;
        g->font = font;
        g->has_titles = has_titles;
        g->ml = 5.5 * font;     /* y tick labels (up to 6 characters) and the rotated y title */
        g->mr = 1.5 * font;     /* half of the last x tick label */
        g->mt = (has_titles ? 2.3 : 0.8) * font;
        g->mb = 3.4 * font;     /* x tick labels and the x title */
        g->gap_x = 3.0 * font;
        g->gap_y = g->rows > 1 ? 3.5 * font : 0.0;   /* room for the next panel title and this panel's last tick label */
        g->pw = (w - g->ml - g->mr - (g->cols - 1) * g->gap_x) / g->cols;
        g->ph = (h - g->mt - g->mb - (g->rows - 1) * g->gap_y) / g->rows;
}

int svg_grid_panel(tld_strbuf* o, const svg_grid* g, int r, int c, const char* title,
                   const svg_axis* x, const svg_axis* y, svg_panel* p)
{
        return svg_grid_panel_ex(o, g, r, c, title, x, y, 0.0, p);
}

int svg_grid_panel_ex(tld_strbuf* o, const svg_grid* g, int r, int c, const char* title,
                      const svg_axis* x, const svg_axis* y, double x_tick_max, svg_panel* p)
{
        double xv[SVG_MAX_TICKS], yv[SVG_MAX_TICKS], xstep = 0.0, ystep = 0.0;
        int xm[SVG_MAX_TICKS], ym[SVG_MAX_TICKS];
        char lab[48];
        double f = g->font;

        if (r < 0 || r >= g->rows || c < 0 || c >= g->cols || g->pw < 20.0 || g->ph < 20.0) {
                ERROR_MSG("svg_plot: panel (%d,%d) does not fit", r, c);
        }
        p->font = f;
        p->w = g->pw;
        p->h = g->ph;
        p->left = g->ml + c * (g->pw + g->gap_x);
        p->top = g->mt + r * (g->ph + g->gap_y);
        p->x = sane(x);
        p->y = sane(y);

        double right = p->left + p->w;
        double bottom = p->top + p->h;
        int nx = svg_ticks(&p->x, (int)(p->w / (4.5 * f)), xv, xm, SVG_MAX_TICKS, &xstep);
        int ny = svg_ticks(&p->y, (int)(p->h / (3.0 * f)), yv, ym, SVG_MAX_TICKS, &ystep);
        if (nx < 0 || ny < 0) {
                ERROR_MSG("svg_plot: ticks failed");
        }
        if (x_tick_max != 0.0 && !p->x.log) {
                /* e.g. an identity axis that runs to 101 so that the 100 % bin is visible */
                while (nx > 0 && xv[nx - 1] > x_tick_max + fabs(xstep) * 1e-9) {
                        nx--;
                }
        }

        /* horizontal grid lines */
        TRY(pf(o, "<path fill=\"none\" stroke=\"" SVG_GRID_COLOR "\" stroke-width=\"1\" d=\""));
        for (int i = 0; i < ny; i++) {
                if (ym[i]) {
                        double yy = py(p, yv[i]);
                        TRY(pf(o, "M%f,%fH%f", p->left, yy, right));
                }
        }
        TRY(pf(o, "\"/>"));

        /* axis lines and tick marks */
        TRY(pf(o, "<path fill=\"none\" stroke=\"" SVG_AXIS_COLOR "\" stroke-width=\"1\" d=\"M%f,%fV%fH%f",
               p->left, p->top, bottom, right));
        for (int i = 0; i < nx; i++) {
                double xx = px(p, xv[i]);
                TRY(pf(o, "M%f,%fv%s", xx, bottom, xm[i] ? "4" : "2.5"));
        }
        for (int i = 0; i < ny; i++) {
                double yy = py(p, yv[i]);
                TRY(pf(o, "M%f,%fh-%s", p->left, yy, ym[i] ? "4" : "2.5"));
        }
        TRY(pf(o, "\"/>"));

        if (title && *title) {
                /* a long title starts further left rather than running off the figure */
                double tx = p->left;
                double tw = 0.62 * f * (double)strlen(title);
                if (tx + tw > g->width - 4.0) {
                        tx = g->width - 4.0 - tw;
                }
                if (tx < 4.0) {
                        tx = 4.0;
                }
                TRY(pf(o, "<text x=\"%f\" y=\"%f\" font-weight=\"bold\">%e</text>", tx, p->top - 0.85 * f, title));
        }
        if (r == g->rows - 1) {
                for (int i = 0; i < nx; i++) {
                        if (xm[i]) {
                                svg_tick_label(&p->x, xstep, xv[i], lab, sizeof(lab));
                                TRY(pf(o, "<text x=\"%f\" y=\"%f\" text-anchor=\"middle\">%e</text>",
                                       px(p, xv[i]), bottom + 7.0 + 0.75 * f, lab));
                        }
                }
                if (p->x.label && *p->x.label) {
                        /* centred under the plot area, moved left or right if it would leave the figure */
                        double half = 0.29 * f * (double)strlen(p->x.label);
                        double cx = p->left + p->w / 2.0;
                        if (cx + half > g->width - 4.0) {
                                cx = g->width - 4.0 - half;
                        }
                        if (cx - half < 4.0) {
                                cx = 4.0 + half;
                        }
                        TRY(pf(o, "<text x=\"%f\" y=\"%f\" text-anchor=\"middle\">%e</text>",
                               cx, bottom + 7.0 + 2.3 * f, p->x.label));
                }
        }
        if (c == 0) {
                for (int i = 0; i < ny; i++) {
                        if (ym[i]) {
                                svg_tick_label(&p->y, ystep, yv[i], lab, sizeof(lab));
                                TRY(pf(o, "<text x=\"%f\" y=\"%f\" text-anchor=\"end\">%e</text>",
                                       p->left - 7.0, py(p, yv[i]) + 0.35 * f, lab));
                        }
                }
                if (p->y.label && *p->y.label) {
                        TRY(pf(o, "<text transform=\"translate(%f %f) rotate(-90)\" text-anchor=\"middle\">%e</text>",
                               0.95 * f, p->top + p->h / 2.0, p->y.label));
                }
        }
        return OK;
ERROR:
        return FAIL;
}

/* ------------------------------------------------------------------ */
/* series                                                               */
/* ------------------------------------------------------------------ */

static int series_open(tld_strbuf* o, const svg_series* s, double width)
{
        TRY(pf(o, "<path fill=\"none\" stroke=\"%e\" stroke-width=\"%f\" stroke-linejoin=\"round\"",
               s->color ? s->color : "#000000", width));
        if (s->dash && *s->dash) {
                TRY(pf(o, " stroke-dasharray=\"%e\"", s->dash));
        }
        TRY(pf(o, " d=\""));
        return OK;
ERROR:
        return FAIL;
}

static int series_close(tld_strbuf* o, const svg_series* s)
{
        if (s->title && *s->title) {
                return pf(o, "\"><title>%e</title></path>", s->title);
        }
        return pf(o, "\"/>");
}

int svg_line(tld_strbuf* o, const svg_panel* p, const svg_series* s)
{
        int valid = 0;
        int run = 0;
        double pend_x = 0.0, pend_y = 0.0;

        if (!s || (s->n > 0 && !s->y)) {
                ERROR_MSG("svg_plot: bad series");
        }
        for (int i = 0; i < s->n; i++) {
                valid += isfinite(s->y[i]) && isfinite(series_x(s, i));
        }
        if (!valid) {
                return OK;
        }
        TRY(series_open(o, s, SVG_LINE_WIDTH));
        for (int i = 0; i <= s->n; i++) {
                int ok = i < s->n && isfinite(s->y[i]) && isfinite(series_x(s, i));
                if (ok) {
                        double xx = px(p, series_x(s, i));
                        double yy = py(p, s->y[i]);
                        if (run == 0) {
                                pend_x = xx;
                                pend_y = yy;
                        } else if (run == 1) {
                                TRY(pf(o, "M%f,%f %f,%f", pend_x, pend_y, xx, yy));
                        } else {
                                TRY(pf(o, " %f,%f", xx, yy));
                        }
                        run++;
                } else {
                        if (run == 1) {
                                /* an isolated value between gaps: a 1.6 px dash so it stays visible */
                                TRY(pf(o, "M%f,%fH%f", pend_x - 0.8, pend_y, pend_x + 0.8));
                        }
                        run = 0;
                }
        }
        return series_close(o, s);
ERROR:
        return FAIL;
}

/* Is bin i of a step series drawn? */
static int step_drawn(const svg_panel* p, const svg_series* s, int i)
{
        return i >= 0 && i < s->n && isfinite(s->y[i]) && s->y[i] > p->y.lo &&
               series_x(s, i + 1) > p->x.lo && series_x(s, i) < p->x.hi;
}

int svg_step(tld_strbuf* o, const svg_panel* p, const svg_series* s)
{
        int any = 0;
        int in_run = 0;
        double base = p->top + p->h;
        double last_y = 0.0;

        if (!s || (s->n > 0 && !s->y)) {
                ERROR_MSG("svg_plot: bad series");
        }
        for (int i = 0; i < s->n && !any; i++) {
                any = step_drawn(p, s, i);
        }
        if (!any) {
                return OK;
        }
        TRY(series_open(o, s, SVG_LINE_WIDTH));
        /*
         * One closed-to-the-baseline outline per run of non-empty bins:
         * up, across, up or down, across, ..., down. Empty bins are skipped.
         * A non-empty bin is at least 1 unit high, and a bin with no drawn neighbour is at least
         * SVG_MIN_BAR wide (centred on the bin), so a lone read on a log axis is still visible.
         */
        for (int i = 0; i < s->n; i++) {
                if (step_drawn(p, s, i)) {
                        double xl = px(p, series_x(s, i));
                        double xr = px(p, series_x(s, i + 1));
                        double yy = py(p, s->y[i]);
                        if (base - yy < 1.0) {
                                yy = base - 1.0;
                        }
                        if (xr - xl < SVG_MIN_BAR && !step_drawn(p, s, i - 1) && !step_drawn(p, s, i + 1)) {
                                double c = (xl + xr) / 2.0;
                                xl = c - SVG_MIN_BAR / 2.0;
                                xr = c + SVG_MIN_BAR / 2.0;
                                if (xl < p->left) {
                                        xr += p->left - xl;
                                        xl = p->left;
                                }
                                if (xr > p->left + p->w) {
                                        xl -= xr - (p->left + p->w);
                                        xr = p->left + p->w;
                                }
                        }
                        if (!in_run) {
                                TRY(pf(o, "M%f,%fV%fH%f", xl, base, yy, xr));
                                in_run = 1;
                        } else if (fabs(yy - last_y) > 0.049) {
                                TRY(pf(o, "V%fH%f", yy, xr));
                        } else {
                                TRY(pf(o, "H%f", xr));
                        }
                        last_y = yy;
                } else if (in_run) {
                        TRY(pf(o, "V%f", base));
                        in_run = 0;
                }
        }
        if (in_run) {
                TRY(pf(o, "V%f", base));
        }
        return series_close(o, s);
ERROR:
        return FAIL;
}

int svg_panel_note(tld_strbuf* o, const svg_panel* p, const char* text)
{
        if (!text || !*text) {
                return OK;
        }
        return pf(o, "<text x=\"%f\" y=\"%f\" text-anchor=\"middle\" fill=\"" SVG_AXIS_COLOR "\">%e</text>",
                  p->left + p->w / 2.0, p->top + p->h / 2.0 + 0.35 * p->font, text);
}

int svg_vline(tld_strbuf* o, const svg_panel* p, double x, const char* label)
{
        double f = p->font;
        double xx = px(p, x);
        int left_side = xx > p->left + 0.6 * p->w;
        double tx = left_side ? xx - 5.0 : xx + 5.0;
        double ty = p->top + 1.1 * f;
        const char* anchor = left_side ? "end" : "start";

        TRY(pf(o, "<path fill=\"none\" stroke=\"#000000\" stroke-width=\"1.25\" stroke-dasharray=\"5 3\" d=\"M%f,%fV%f\">",
               xx, p->top, p->top + p->h));
        TRY(pf(o, "<title>%e</title></path>", label));
        if (label && *label) {
                /* white halo drawn first so the label stays readable over lines */
                TRY(pf(o, "<text x=\"%f\" y=\"%f\" text-anchor=\"%s\" fill=\"#fff\" stroke=\"#fff\" stroke-width=\"4\" "
                       "stroke-linejoin=\"round\">%e</text>", tx, ty, anchor, label));
                TRY(pf(o, "<text x=\"%f\" y=\"%f\" text-anchor=\"%s\">%e</text>", tx, ty, anchor, label));
        }
        return OK;
ERROR:
        return FAIL;
}

/* ------------------------------------------------------------------ */
/* swatches                                                             */
/* ------------------------------------------------------------------ */

static int sample(tld_strbuf* o, double x0, double x1, double y, const char* color, const char* dash)
{
        TRY(pf(o, "<path fill=\"none\" stroke=\"%e\" stroke-width=\"2.5\"", color));
        if (dash && *dash) {
                TRY(pf(o, " stroke-dasharray=\"%e\"", dash));
        }
        TRY(pf(o, " d=\"M%f,%fH%f\"/>", x0, y, x1));
        return OK;
ERROR:
        return FAIL;
}

int svg_swatch(tld_strbuf* o, const char* id, const char* color, const char* dash,
               const char* title, const char* desc)
{
        TRY(begin_ex(o, id, 52.0, 14.0, SVG_DEFAULT_FONT, title, desc, "vertical-align:middle"));
        TRY(sample(o, 2.0, 50.0, 7.0, color, dash));
        return svg_end(o);
ERROR:
        return FAIL;
}

int svg_swatch_row(tld_strbuf* o, const char* id, const svg_swatch_item* items, int n,
                   double font, const char* title, const char* desc)
{
        double h = 1.8 * font;
        double x = 2.0;

        for (int i = 0; i < n; i++) {
                x += 36.0 + 6.0 + 0.62 * font * (double)strlen(items[i].label ? items[i].label : "") + 14.0;
        }
        TRY(begin_ex(o, id, x, h, font, title, desc, "max-width:100%;height:auto;vertical-align:top"));
        x = 2.0;
        for (int i = 0; i < n; i++) {
                const char* label = items[i].label ? items[i].label : "";
                TRY(sample(o, x, x + 36.0, h / 2.0, items[i].color, items[i].dash));
                TRY(pf(o, "<text x=\"%f\" y=\"%f\">%e</text>", x + 42.0, h / 2.0 + 0.35 * font, label));
                x += 42.0 + 0.62 * font * (double)strlen(label) + 14.0;
        }
        return svg_end(o);
ERROR:
        return FAIL;
}
