/*
 * Unit test for the SVG primitives (WP1.6).
 *
 * Renders every primitive into memory and checks: tick logic and labels,
 * exact path data on hand-computable cases, the encoding rules (one <path> and
 * one <title> per series, no per-point elements), deterministic bytes, size
 * bounds, and well-formedness of every SVG (see svg_check).
 *
 * With one argument the specimen sheet is also written there as an HTML page,
 * for looking at the figures in a browser:  svg_plot_test /some/dir/specimen.html
 * (ctest runs it without arguments and writes nothing).
 */

#include "svg_plot.h"

#include <ctype.h>
#include <locale.h>
#include <math.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond, ...)                                                \
        do {                                                            \
                if (!(cond)) {                                          \
                        failures++;                                     \
                        fprintf(stderr, "FAIL %s:%d: %s\n  ", __FILE__, __LINE__, #cond); \
                        fprintf(stderr, __VA_ARGS__);                   \
                        fprintf(stderr, "\n");                          \
                }                                                       \
        } while (0)

#include "svg_check.h"

/* ------------------------------------------------------------------ */
/* helpers                                                              */
/* ------------------------------------------------------------------ */

static int count_of(const char* hay, const char* needle)
{
        int n = 0;
        size_t l = strlen(needle);
        for (const char* p = hay; (p = strstr(p, needle)); p += l) {
                n++;
        }
        return n;
}

/* Copy the d="..." of the path whose series-level <title> is title. */
static int path_d_of(const char* doc, const char* title, char* d, size_t dn)
{
        char t[128];
        snprintf(t, sizeof(t), "<title>%s</title>", title);
        const char* at = strstr(doc, t);
        if (!at) {
                return 0;
        }
        const char* start = NULL;
        for (const char* p = doc; p < at; p++) {
                if (strncmp(p, " d=\"", 4) == 0) {
                        start = p + 4;
                }
        }
        if (!start) {
                return 0;
        }
        size_t n = (size_t)(strchr(start, '"') - start);
        if (n >= dn) {
                n = dn - 1;
        }
        memcpy(d, start, n);
        d[n] = 0;
        return 1;
}

static uint32_t lcg_state = 12345;
static double rnd(void)
{
        lcg_state = lcg_state * 1664525u + 1013904223u;
        return (lcg_state >> 8) / 16777216.0;
}

/* ------------------------------------------------------------------ */
/* tests                                                                */
/* ------------------------------------------------------------------ */

static void test_ticks(void)
{
        double v[SVG_MAX_TICKS], step = 0;
        int major[SVG_MAX_TICKS];
        char lab[48];

        CHECK(svg_nice_step(100, 5) == 20.0, "%g", svg_nice_step(100, 5));
        CHECK(svg_nice_step(1, 5) == 0.2, "%g", svg_nice_step(1, 5));
        CHECK(svg_nice_step(1000, 4) == 200.0, "%g", svg_nice_step(1000, 4));
        CHECK(svg_nice_step(0, 5) == 1.0, "range 0");
        CHECK(svg_nice_step(NAN, 5) == 1.0, "range nan");
        CHECK(svg_nice_ceil(173, 5) == 200.0, "%g", svg_nice_ceil(173, 5));
        CHECK(svg_nice_ceil(100, 5) == 100.0, "%g", svg_nice_ceil(100, 5));
        CHECK(svg_nice_ceil(0, 5) == 1.0, "max 0");

        svg_axis lin = {0, 150, 0, NULL};
        int n = svg_ticks(&lin, 6, v, major, SVG_MAX_TICKS, &step);
        CHECK(n == 8 && step == 20.0 && v[0] == 0.0 && v[7] == 140.0, "n=%d step=%g", n, step);

        svg_axis unit = {0, 1, 0, NULL};
        n = svg_ticks(&unit, 5, v, major, SVG_MAX_TICKS, &step);
        CHECK(n == 6 && v[5] == 1.0 && step == 0.2, "n=%d step=%g", n, step);
        svg_tick_label(&unit, step, v[1], lab, sizeof(lab));
        CHECK(strcmp(lab, "0.2") == 0, "%s", lab);

        /* many decades: decades labelled, 2 and 5 as unlabelled ticks */
        svg_axis lg = {100, 1e6, 1, NULL};
        n = svg_ticks(&lg, 8, v, major, SVG_MAX_TICKS, &step);
        int nmaj = 0;
        for (int i = 0; i < n; i++) {
                nmaj += major[i];
        }
        CHECK(n == 13 && nmaj == 5 && step == 0.0, "n=%d major=%d", n, nmaj);
        const char* want[] = {"100", "1 k", "10 k", "100 k", "1 M"};
        int w = 0;
        for (int i = 0; i < n; i++) {
                if (major[i]) {
                        svg_tick_label(&lg, step, v[i], lab, sizeof(lab));
                        CHECK(strcmp(lab, want[w]) == 0, "log label %d: %s", w, lab);
                        w++;
                }
        }
        /* few decades: 1, 2, 5 all labelled */
        svg_axis lg2 = {100, 3000, 1, NULL};
        n = svg_ticks(&lg2, 12, v, major, SVG_MAX_TICKS, &step);
        CHECK(n == 5 && v[0] == 100 && v[1] == 200 && v[2] == 500 && v[3] == 1000 && v[4] == 2000 &&
              major[1] == 1, "n=%d", n);
        /* a short log axis (a lone bin: 1.4 decades) still has at least three labelled ticks */
        svg_axis lg1 = {200000, 5e6, 1, NULL};
        n = svg_ticks(&lg1, 4, v, major, SVG_MAX_TICKS, &step);
        nmaj = 0;
        for (int i = 0; i < n; i++) {
                nmaj += major[i];
        }
        CHECK(n == 5 && nmaj == 5 && v[0] == 2e5 && v[2] == 1e6 && v[4] == 5e6, "short log axis: n=%d major=%d", n, nmaj);
        /* too many decades for the room: every third decade */
        svg_axis lg3 = {1, 1e12, 1, NULL};
        n = svg_ticks(&lg3, 5, v, major, SVG_MAX_TICKS, &step);
        CHECK(n == 5 && v[1] == 1e3 && v[4] == 1e12, "n=%d", n);
        /* degenerate input is repaired, never crashes */
        svg_axis bad = {5, 5, 1, NULL};
        CHECK(svg_ticks(&bad, 5, v, major, SVG_MAX_TICKS, &step) > 0, "degenerate log");
        CHECK(svg_ticks(NULL, 5, v, major, SVG_MAX_TICKS, &step) == -1, "null axis");

        /* linear label units */
        svg_axis a4 = {0, 40000, 0, NULL};
        svg_tick_label(&a4, 10000, 20000, lab, sizeof(lab));
        CHECK(strcmp(lab, "20 k") == 0, "%s", lab);
        svg_tick_label(&a4, 10000, 0, lab, sizeof(lab));
        CHECK(strcmp(lab, "0") == 0, "zero has no unit: %s", lab);
        svg_axis a6 = {0, 1.8e6, 0, NULL};
        svg_tick_label(&a6, 5e5, 1.5e6, lab, sizeof(lab));
        CHECK(strcmp(lab, "1.5 M") == 0, "%s", lab);
        svg_axis a1 = {0, 1000, 0, NULL};
        svg_tick_label(&a1, 200, 400, lab, sizeof(lab));
        CHECK(strcmp(lab, "400") == 0, "%s", lab);
        svg_tick_label(&a1, 200, -1e-15, lab, sizeof(lab));
        CHECK(strcmp(lab, "0") == 0, "negative zero: %s", lab);
}

static const svg_axis AX_X5 = {0, 5, 0, NULL};
static const svg_axis AX_Y10 = {0, 10, 0, NULL};

static svg_panel unit_panel(void)
{
        svg_panel p = {0, 0, 100, 100, 12, AX_X5, AX_Y10};
        return p;
}

static void test_series_paths(void)
{
        tld_strbuf* o = NULL;
        char d[512];
        svg_panel p = unit_panel();
        tld_strbuf_alloc(&o, 1024);

        /* line with gaps: runs [1,2], [4,5] and an isolated 7 */
        double y[] = {1, 2, NAN, 4, NAN, 7};
        svg_series gap = {"gap series", "#0072B2", "6 3", NULL, 0, 1, y, 6};
        CHECK(svg_line(o, &p, &gap) == OK, "line");
        CHECK(path_d_of(TLD_STR(o), "gap series", d, sizeof(d)), "no path");
        CHECK(strcmp(d, "M0,90 20,80M59.2,60H60.8M99.2,30H100.8") == 0, "%s", d);
        CHECK(strstr(TLD_STR(o), "stroke-dasharray=\"6 3\"") != NULL, "dash");
        CHECK(count_of(TLD_STR(o), "<path") == 1 && count_of(TLD_STR(o), "<title>") == 1, "one path, one title");
        tld_strbuf_clear(o);

        /* exact coordinates: x 0..5 over 100 px, y 0..10 over 100 px */
        double yl[] = {0, 5, 10};
        svg_series ln = {"ln", "#000000", NULL, NULL, 0, 1, yl, 3};
        CHECK(svg_line(o, &p, &ln) == OK, "line");
        CHECK(path_d_of(TLD_STR(o), "ln", d, sizeof(d)) && strcmp(d, "M0,100 20,50 40,0") == 0, "%s", d);
        CHECK(strstr(TLD_STR(o), "stroke-dasharray") == NULL, "solid has no dasharray");
        tld_strbuf_clear(o);

        /* explicit x */
        double xe[] = {0.5, 1.5, 4.5};
        svg_series lx = {"lx", "#000000", "", xe, 0, 0, yl, 3};
        CHECK(svg_line(o, &p, &lx) == OK, "line");
        CHECK(path_d_of(TLD_STR(o), "lx", d, sizeof(d)) && strcmp(d, "M10,100 30,50 90,0") == 0, "%s", d);
        tld_strbuf_clear(o);

        /* all gaps: nothing drawn */
        double ng[] = {NAN, NAN};
        svg_series none = {"none", "#000000", NULL, NULL, 0, 1, ng, 2};
        CHECK(svg_line(o, &p, &none) == OK && TLD_STRLEN(o) == 0, "all-gap series draws nothing");

        /* out-of-range values are clamped into the plotting area, never leak out of it */
        double big[] = {-5, 100, INFINITY};
        svg_series cl = {"cl", "#000000", NULL, NULL, 0, 1, big, 3};
        CHECK(svg_line(o, &p, &cl) == OK, "clamp");
        CHECK(path_d_of(TLD_STR(o), "cl", d, sizeof(d)) && strcmp(d, "M0,100 20,0") == 0, "%s", d);
        tld_strbuf_clear(o);

        /* step histogram: two runs, equal neighbours share one level */
        double edges[] = {0, 1, 2, 3, 4, 5};
        double hy[] = {0, 3, 3, 0, 5};
        svg_series st = {"st", "#D55E00", "2 3", edges, 0, 0, hy, 5};
        CHECK(svg_step(o, &p, &st) == OK, "step");
        CHECK(path_d_of(TLD_STR(o), "st", d, sizeof(d)) &&
              strcmp(d, "M20,100V70H40H60V100M80,100V50H100V100") == 0, "%s", d);
        tld_strbuf_clear(o);

        /* step with uniform edges, a NaN bin and a bin that touches the end */
        double hy2[] = {2, NAN, 4, 4, 6};
        svg_series st2 = {"st2", "#000000", NULL, NULL, 0, 1, hy2, 5};
        CHECK(svg_step(o, &p, &st2) == OK, "step2");
        CHECK(path_d_of(TLD_STR(o), "st2", d, sizeof(d)) &&
              strcmp(d, "M0,100V80H20V100M40,100V60H60H80V40H100V100") == 0, "%s", d);
        tld_strbuf_clear(o);

        /* a lone, very narrow bin (0.1 px wide, 0.5 px high) is drawn at least 6 px wide and 1 px high,
         * so a single read on a log axis is not a hairline; bins with neighbours keep their own edges */
        {
                static double lone[1000];
                double xl = 0, xr = 0, top = 0, base = 0;
                memset(lone, 0, sizeof(lone));
                lone[500] = 0.05;
                svg_series ln1 = {"lone", "#000000", NULL, NULL, 0, 0.005, lone, 1000};
                CHECK(svg_step(o, &p, &ln1) == OK, "lone step");
                CHECK(path_d_of(TLD_STR(o), "lone", d, sizeof(d)) &&
                      sscanf(d, "M%lf,%lfV%lfH%lfV", &xl, &base, &top, &xr) == 4, "lone path: %s", d);
                CHECK(xr - xl >= 5.99, "lone bin is %g px wide: %s", xr - xl, d);
                CHECK(base - top >= 0.99, "lone bin is %g px high: %s", base - top, d);
                CHECK((xl + xr) / 2.0 > 49.9 && (xl + xr) / 2.0 < 50.2, "lone bin stays centred: %s", d);
                tld_strbuf_clear(o);
                lone[501] = 0.05;       /* a neighbour: the bins form one run and keep their exact width */
                CHECK(svg_step(o, &p, &ln1) == OK, "pair step");
                CHECK(path_d_of(TLD_STR(o), "lone", d, sizeof(d)) &&
                      sscanf(d, "M%lf,%lfV%lfH%lfV", &xl, &base, &top, &xr) == 4 && xr - xl < 0.5, "run keeps exact width: %s", d);
                tld_strbuf_clear(o);
        }

        /* log x axis (plot 1): edge 0 clamps to the axis start, decades are equally spaced */
        svg_panel lp = {0, 0, 300, 100, 12, {1, 1000, 1, NULL}, AX_Y10};
        double lx2[] = {1, 10, 100, 1000};
        double ly[] = {0, 5, 10, 0};
        svg_series ll = {"ll", "#000000", NULL, lx2, 0, 0, ly, 4};
        CHECK(svg_line(o, &lp, &ll) == OK, "log line");
        CHECK(path_d_of(TLD_STR(o), "ll", d, sizeof(d)) && strcmp(d, "M0,100 100,50 200,0 300,100") == 0, "%s", d);
        tld_strbuf_clear(o);
        double e0[] = {0, 10, 100, 1000};
        double h0[] = {4, 4, 8};
        svg_series ls = {"ls", "#000000", NULL, e0, 0, 0, h0, 3};
        CHECK(svg_step(o, &lp, &ls) == OK, "log step");
        CHECK(path_d_of(TLD_STR(o), "ls", d, sizeof(d)) && strncmp(d, "M0,100V60H100H", 14) == 0, "%s", d);
        tld_strbuf_clear(o);

        /* a series without a title still gets a path */
        svg_series nt = {NULL, "#000000", NULL, NULL, 0, 1, yl, 3};
        CHECK(svg_line(o, &p, &nt) == OK && count_of(TLD_STR(o), "<title>") == 0 &&
              count_of(TLD_STR(o), "<path") == 1, "untitled series");
        tld_strbuf_clear(o);

        /* vertical marker */
        CHECK(svg_vline(o, &p, 2.5, "N50 12 k") == OK, "vline");
        CHECK(strstr(TLD_STR(o), "d=\"M50,0V100\"") != NULL, "%s", TLD_STR(o));
        CHECK(strstr(TLD_STR(o), "stroke-dasharray=\"5 3\"") != NULL, "marker dash");
        CHECK(count_of(TLD_STR(o), ">N50 12 k</text>") == 2, "label (halo + text)");
        CHECK(strstr(TLD_STR(o), "text-anchor=\"start\"") != NULL, "left of centre: label on the right");
        tld_strbuf_clear(o);
        CHECK(svg_vline(o, &p, 4.5, "median 98.2%") == OK &&
              strstr(TLD_STR(o), "text-anchor=\"end\"") != NULL, "right of centre: label on the left");

        tld_strbuf_free(o);
}

static void test_escaping_and_errors(void)
{
        tld_strbuf* o = NULL;
        char err[200];
        tld_strbuf_alloc(&o, 1024);

        CHECK(svg_begin(o, "e1", 100, 50, 12, "A & B <x> \"q\" 'r'", "line\nbreak") == OK, "begin");
        CHECK(svg_end(o) == OK, "end");
        CHECK(strstr(TLD_STR(o), "A &amp; B &lt;x&gt; &quot;q&quot; &#39;r&#39;") != NULL, "%s", TLD_STR(o));
        CHECK(strstr(TLD_STR(o), "line break") != NULL, "control characters become spaces");
        CHECK(svg_check(TLD_STR(o), err, sizeof(err)), "%s", err);
        tld_strbuf_clear(o);

        /* bad parameters fail cleanly (error messages on stderr are expected here) */
        CHECK(svg_begin(o, "e2", 100, 50, 10, "t", "d") == FAIL, "font below 12 px is refused");
        CHECK(svg_begin(o, "", 100, 50, 12, "t", "d") == FAIL, "empty id is refused");
        CHECK(svg_begin(o, "e3", 0, 50, 12, "t", "d") == FAIL, "zero width is refused");
        tld_strbuf_clear(o);

        svg_grid g;
        svg_grid_init(&g, 1, 1, 30, 30, 12, 0);
        svg_panel p;
        svg_axis ax = {0, 1, 0, NULL};
        CHECK(svg_grid_panel(o, &g, 0, 0, NULL, &ax, &ax, &p) == FAIL, "a panel that does not fit is refused");
        svg_grid_init(&g, 1, 1, 400, 200, 12, 0);
        CHECK(svg_grid_panel(o, &g, 1, 0, NULL, &ax, &ax, &p) == FAIL, "row out of range");
        tld_strbuf_free(o);
}

static void test_checker(void)
{
        char err[200];
        const char* ok = "<svg role=\"img\"><title>t</title><desc>d</desc><path d=\"M0,0\"/></svg>\n";
        CHECK(svg_check(ok, err, sizeof(err)), "%s", err);

        const char* bad[] = {
            "<svg role=\"img\"><title>t</title><desc>d</desc><g></svg>",              /* unbalanced */
            "<svg role=\"img\"><title>t</title><desc>d</desc></g></svg>",             /* stray close */
            "<svg role=\"img\"><title>t</title><desc>d</desc><path d=M0,0/></svg>",   /* unquoted */
            "<svg role=\"img\"><title>t</title><desc>d</desc><path d=\"M0,0/></svg>", /* unterminated */
            "<svg role=\"img\"><title>t</title><desc>d</desc><g a=\"1\" a=\"2\"/></svg>", /* duplicate */
            "<svg role=\"img\"><title>A & B</title><desc>d</desc></svg>",             /* raw ampersand */
            "<svg role=\"img\"><title>t</title><desc>d</desc><path d=\"M0,nan\"/></svg>", /* non-finite */
            "<svg role=\"img\"><title>t</title><desc>d</desc><path d=\"M0,inf\"/></svg>",
            "<svg><title>t</title><desc>d</desc></svg>",                              /* no role */
            "<svg role=\"img\"><desc>d</desc></svg>",                                 /* no title */
            "<svg role=\"img\"><title>t</title></svg>",                               /* no desc */
            "<svg role=\"img\"><title>t</title><desc>d</desc><circle/></svg>",        /* per-point element */
            "<svg role=\"img\"><title>t</title><desc>d</desc><script>x</script></svg>",
            "<svg role=\"img\"><title id=\"a\">t</title><desc id=\"a\">d</desc></svg>", /* duplicate id */
            "<svg role=\"img\"><title>t</title><desc>d</desc>",                       /* not closed */
            "<svg role=\"img\"><title>t</title><desc>d</desc><!-- x --></svg>",       /* comment */
            "<g></g>",                                                                /* outside svg */
            "",
        };
        for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); i++) {
                CHECK(!svg_check(bad[i], err, sizeof(err)), "checker accepted: %s", bad[i]);
        }
}

/* ------------------------------------------------------------------ */
/* specimen sheet: the figures the report needs, rendered with fake data */
/* ------------------------------------------------------------------ */

static const char* MAPQ_COLOR[4] = {"#0072B2", "#B86E00", "#D55E00", "#777777"};
static const char* MAPQ_DASH[4] = {NULL, "6 3", "2 3", "8 3 2 3"};
static const char* MAPQ_NAME[4] = {"MAPQ >= 30", "MAPQ 10-29", "MAPQ < 10", "Unmapped"};

#define NBIN 1000

static int fig_length(tld_strbuf* o)
{
        enum { N = 60 };
        static double edges[N + 1], y[4][N];
        svg_series s[4];
        svg_grid g;
        svg_panel p;
        double top = 0;

        for (int i = 0; i <= N; i++) {
                edges[i] = 100.0 * pow(10.0, i / 20.0);
        }
        for (int c = 0; c < 4; c++) {
                for (int i = 0; i < N; i++) {
                        double z = (i - 24.0 - 3.0 * c) / 9.0;
                        y[c][i] = (c == 0 ? 90000.0 : 30000.0 / (c + 1)) * exp(-z * z) * (0.9 + 0.2 * rnd());
                        top = fmax(top, y[c][i]);
                }
                s[c] = (svg_series){MAPQ_NAME[c], MAPQ_COLOR[c], MAPQ_DASH[c], edges, 0, 0, y[c], N};
        }
        svg_axis xa = {100, 1e5, 1, "read length (bp)"};
        svg_axis ya = {0, svg_nice_ceil(top, 5), 0, "reads per bin"};
        svg_grid_init(&g, 1, 1, 640, 300, 12, 1);
        RUN(svg_begin(o, "p1", 640, 300, 12, "How long are the reads?",
                      "Reads per length bin on a log scale, one line per mapping-quality class; N50 is marked."));
        RUN(svg_grid_panel(o, &g, 0, 0, NULL, &xa, &ya, &p));
        for (int c = 3; c >= 0; c--) {
                RUN(svg_step(o, &p, &s[c]));
        }
        RUN(svg_vline(o, &p, 11800.0, "N50 11.8 k"));
        return svg_end(o);
ERROR:
        return FAIL;
}

static int fig_quality(tld_strbuf* o, int paired)
{
        static double y[2][3][NBIN];
        svg_series s;
        svg_grid g;
        svg_panel p;
        int cols = paired ? 2 : 1;

        svg_grid_init(&g, 1, cols, paired ? 640 : 480, 280, 12, 1);
        RUN(svg_begin(o, paired ? "p3b" : "p3", g.width, g.height, 12, "Where along the read does quality drop?",
                      "Mean base quality by position, one line per mapping-quality class."));
        svg_axis xa = {1, 150, 0, "position in read (bp), aligned part"};
        svg_axis ya = {0, 40, 0, "mean Q"};
        for (int m = 0; m < cols; m++) {
                RUN(svg_grid_panel(o, &g, 0, m, m ? "Read 2" : (paired ? "Read 1" : "All reads"), &xa, &ya, &p));
                for (int c = 0; c < 3; c++) {
                        for (int i = 0; i < 150; i++) {
                                y[m][c][i] = 37.0 - 0.02 * i * (m + 1) - 2.0 * c + 1.5 * rnd();
                                if (c == 2 && i % 40 > 30) {
                                        y[m][c][i] = NAN;
                                }
                        }
                        s = (svg_series){MAPQ_NAME[c], MAPQ_COLOR[c], MAPQ_DASH[c], NULL, 1, 1, y[m][c], 150};
                        RUN(svg_line(o, &p, &s));
                }
        }
        return svg_end(o);
ERROR:
        return FAIL;
}

static int fig_composition(tld_strbuf* o)
{
        static const svg_swatch_item base[5] = {
            {"A", "#009E73", NULL}, {"C", "#CC79A7", NULL}, {"G", "#000000", NULL},
            {"T", "#8C5A00", NULL}, {"N", "#2E86C1", "1 2"}};
        static double y[2][5][NBIN];
        svg_grid g;
        svg_panel p;

        RUN(svg_swatch_row(o, "p4k", base, 5, 12, "Base colours", "Line colours for A, C, G, T and N."));
        svg_grid_init(&g, 1, 2, 640, 280, 12, 1);
        RUN(svg_begin(o, "p4", 640, 280, 12, "Is base composition even along the read?",
                      "Percentage of A, C, G, T and N by position, for read 1 and read 2."));
        svg_axis xa = {1, 150, 0, "position in read (bp), aligned part"};
        svg_axis ya = {0, 100, 0, "% of bases"};
        for (int m = 0; m < 2; m++) {
                RUN(svg_grid_panel(o, &g, 0, m, m ? "Read 2" : "Read 1", &xa, &ya, &p));
                for (int b = 0; b < 5; b++) {
                        for (int i = 0; i < 150; i++) {
                                y[m][b][i] = b == 4 ? 0.3 * rnd() : 25.0 + 2.0 * (rnd() - 0.5) + (i < 8 ? 6.0 * (b % 2 ? 1 : -1) : 0);
                        }
                        svg_series s = {base[b].label, base[b].color, base[b].dash, NULL, 1, 1, y[m][b], 150};
                        RUN(svg_line(o, &p, &s));
                }
        }
        return svg_end(o);
ERROR:
        return FAIL;
}

static int fig_errors(tld_strbuf* o)
{
        static const char* row_label[3] = {"Mismatches (%)", "Insertions (%)", "Deletions (%)"};
        static const double row_max[3] = {2.0, 0.5, 0.5};
        static double y[3][2][3][NBIN];
        svg_grid g;
        svg_panel p;

        svg_grid_init(&g, 3, 2, 640, 520, 12, 1);
        RUN(svg_begin(o, "p5", 640, 520, 12, "Where along the read do errors occur?",
                      "Mismatch, insertion and deletion rates per aligned base by position, for read 1 and read 2."));
        for (int r = 0; r < 3; r++) {
                svg_axis xa = {1, 150, 0, "position in read (bp), aligned part"};
                svg_axis ya = {0, row_max[r], 0, row_label[r]};
                for (int m = 0; m < 2; m++) {
                        RUN(svg_grid_panel(o, &g, r, m, r == 0 ? (m ? "Read 2" : "Read 1") : NULL, &xa, &ya, &p));
                        for (int c = 0; c < 3; c++) {
                                for (int i = 0; i < 150; i++) {
                                        y[r][m][c][i] = row_max[r] * (0.2 + 0.2 * c + 0.3 * rnd() + (i > 140 ? 0.3 : 0.0));
                                }
                                svg_series s = {MAPQ_NAME[c], MAPQ_COLOR[c], MAPQ_DASH[c], NULL, 1, 1, y[r][m][c], 150};
                                RUN(svg_line(o, &p, &s));
                        }
                }
        }
        return svg_end(o);
ERROR:
        return FAIL;
}

static int fig_identity(tld_strbuf* o)
{
        enum { N = 60 };
        static double y[3][N];
        svg_series s;
        svg_grid g;
        svg_panel p;
        double top = 0;

        svg_grid_init(&g, 1, 1, 480, 280, 12, 0);
        RUN(svg_begin(o, "p6", 480, 280, 12, "How accurate are the aligned reads?",
                      "Reads by identity to the reference, one line per mapping-quality class; the median is marked."));
        for (int c = 0; c < 3; c++) {
                for (int i = 0; i < N; i++) {
                        double z = (i - 52.0 + 3 * c) / 6.0;
                        y[c][i] = (c ? 1500.0 : 9000.0) * exp(-z * z);
                        top = fmax(top, y[c][i]);
                }
        }
        svg_axis xa = {88, 100, 0, "identity (%)"};
        svg_axis ya = {0, svg_nice_ceil(top, 5), 0, "reads"};
        RUN(svg_grid_panel(o, &g, 0, 0, NULL, &xa, &ya, &p));
        for (int c = 0; c < 3; c++) {
                s = (svg_series){MAPQ_NAME[c], MAPQ_COLOR[c], MAPQ_DASH[c], NULL, 88, 0.2, y[c], N};
                RUN(svg_step(o, &p, &s));
        }
        RUN(svg_vline(o, &p, 97.4, "median 97.4%"));
        return svg_end(o);
ERROR:
        return FAIL;
}

static int render_all(tld_strbuf* o)
{
        lcg_state = 12345;
        RUN(fig_length(o));
        RUN(fig_quality(o, 0));
        RUN(fig_quality(o, 1));
        RUN(fig_composition(o));
        RUN(fig_errors(o));
        RUN(fig_identity(o));
        for (int c = 0; c < 4; c++) {
                char id[16];
                snprintf(id, sizeof(id), "sw%d", c);
                RUN(svg_swatch(o, id, MAPQ_COLOR[c], MAPQ_DASH[c], MAPQ_NAME[c], "Line sample"));
        }
        return OK;
ERROR:
        return FAIL;
}

static void test_grid(void)
{
        tld_strbuf* o = NULL;
        svg_grid g;
        svg_panel p00, p21, p10;
        svg_axis xa = {1, 150, 0, "XLAB"}, ya = {0, 1, 0, "YLAB"};
        tld_strbuf_alloc(&o, 1024);

        svg_grid_init(&g, 3, 2, 640, 520, 12, 1);
        CHECK(g.pw > 100 && g.ph > 60, "plotting area %g x %g", g.pw, g.ph);
        svg_begin(o, "g1", 640, 520, 12, "grid", "grid");
        CHECK(svg_grid_panel(o, &g, 0, 0, "Read 1", &xa, &ya, &p00) == OK, "p00");
        CHECK(svg_grid_panel(o, &g, 1, 0, NULL, &xa, &ya, &p10) == OK, "p10");
        CHECK(svg_grid_panel(o, &g, 2, 1, NULL, &xa, &ya, &p21) == OK, "p21");
        svg_end(o);

        CHECK(p00.w == p21.w && p00.h == p21.h, "all plotting areas have the same size");
        CHECK(p21.left > p00.left + p00.w && p10.top > p00.top + p00.h, "panels do not overlap");
        CHECK(p21.top + p21.h + g.mb <= g.height + 1e-9 && p21.left + p21.w + g.mr <= g.width + 1e-9, "fits in the figure");
        CHECK(count_of(TLD_STR(o), ">XLAB<") == 1, "x title only under the bottom row's first call (%d)", count_of(TLD_STR(o), ">XLAB<"));
        CHECK(count_of(TLD_STR(o), ">YLAB<") == 2, "y title in column 0 only (%d)", count_of(TLD_STR(o), ">YLAB<"));
        CHECK(count_of(TLD_STR(o), ">Read 1<") == 1, "panel title");
        /* a note inside a panel: centred in the plotting area, dark enough to read, nothing for NULL or "" */
        tld_strbuf_clear(o);
        CHECK(svg_panel_note(o, &p10, "none detected") == OK && svg_panel_note(o, &p10, "") == OK &&
              svg_panel_note(o, &p10, NULL) == OK, "panel note");
        CHECK(count_of(TLD_STR(o), "<text") == 1 && count_of(TLD_STR(o), ">none detected<") == 1 &&
              strstr(TLD_STR(o), "fill=\"#595959\"") != NULL && strstr(TLD_STR(o), "text-anchor=\"middle\"") != NULL,
              "%s", TLD_STR(o));
        tld_strbuf_free(o);
}

static void test_tick_max(void)
{
        /* an identity axis that runs to 101 so the 100 % bin shows, without a "101" label */
        tld_strbuf* o = NULL;
        svg_grid g;
        svg_panel p;
        svg_axis xa = {98, 101, 0, NULL}, ya = {0, 1, 0, NULL};
        char err[200];
        tld_strbuf_alloc(&o, 1024);
        svg_grid_init(&g, 1, 1, 420, 300, 14, 0);
        svg_begin(o, "tm", 420, 300, 14, "t", "d");
        CHECK(svg_grid_panel(o, &g, 0, 0, NULL, &xa, &ya, &p) == OK, "panel");
        CHECK(count_of(TLD_STR(o), ">101.0<") == 1 || count_of(TLD_STR(o), ">101<") == 1, "without a cap the 101 tick is labelled");
        tld_strbuf_clear(o);
        svg_begin(o, "tm2", 420, 300, 14, "t", "d");
        CHECK(svg_grid_panel_ex(o, &g, 0, 0, NULL, &xa, &ya, 100.0, &p) == OK, "panel_ex");
        svg_end(o);
        CHECK(strstr(TLD_STR(o), ">101") == NULL && strstr(TLD_STR(o), ">100") != NULL, "the 101 tick is dropped, 100 stays");
        CHECK(p.x.hi == 101.0, "the axis itself still runs to 101");
        CHECK(svg_check(TLD_STR(o), err, sizeof(err)), "%s", err);

        /* rows are far enough apart for the next panel title (1 row of titles + the last tick label) */
        svg_grid_init(&g, 3, 1, 420, 400, 14, 1);
        CHECK(g.gap_y >= 3.0 * 14, "row gap %g", g.gap_y);
        tld_strbuf_free(o);
}

static void test_specimens(const char* dump)
{
        tld_strbuf *a = NULL, *b = NULL;
        char err[300];
        tld_strbuf_alloc(&a, 1 << 16);
        tld_strbuf_alloc(&b, 1 << 16);

        CHECK(render_all(a) == OK, "render_all");
        CHECK(render_all(b) == OK, "render_all (second)");
        CHECK(TLD_STRLEN(a) == TLD_STRLEN(b) && strcmp(TLD_STR(a), TLD_STR(b)) == 0, "output is not deterministic");
        CHECK(svg_check(TLD_STR(a), err, sizeof(err)), "%s", err);

        const char* doc = TLD_STR(a);
        CHECK(count_of(doc, "<svg ") == 11, "svg count %d", count_of(doc, "<svg "));
        CHECK(count_of(doc, "role=\"img\"") == 11, "role=img on every svg");
        CHECK(count_of(doc, "<desc id=") == 11 && count_of(doc, "<title id=") == 11, "title and desc on every svg");
        CHECK(count_of(doc, "aria-labelledby=") == 11, "aria-labelledby on every svg");
        CHECK(strstr(doc, "<script") == NULL && count_of(doc, "http") == count_of(doc, "http://www.w3.org/2000/svg"),
              "no script and no external references (only the SVG namespace)");
        CHECK(strstr(doc, "font-size=\"11") == NULL && strstr(doc, "font-size=\"10") == NULL, "no text below 12 px");

        /* every series is one titled <path>; there are no per-point or loose line elements */
        CHECK(count_of(doc, "<title>MAPQ &gt;= 30</title>") >= 6, "series titles present");
        CHECK(count_of(doc, "<circle") == 0 && count_of(doc, "<line") == 0, "no per-point or loose line elements");

        /* size: no per-point elements, so bounded (the report must stay under 1 MB) */
        size_t total = (size_t)TLD_STRLEN(a);
        CHECK(total < 200000, "specimen sheet is %zu bytes", total);
        fprintf(stderr, "specimen sheet: 11 svgs, %zu bytes\n", total);

        if (dump) {
                FILE* f = fopen(dump, "w");
                if (f) {
                        fprintf(f, "<!doctype html><meta charset=\"utf-8\"><title>svg_plot specimens</title>"
                                "<body style=\"font-family:sans-serif;margin:16px\"><p>Specimens (fake data).</p>%s", doc);
                        fclose(f);
                }
        }
        tld_strbuf_free(a);
        tld_strbuf_free(b);
}

static void test_locale_and_size(void)
{
        /* numbers must not depend on the locale (a decimal comma would corrupt path data) */
        const char* loc = setlocale(LC_ALL, "de_DE.UTF-8");
        if (!loc) {
                loc = setlocale(LC_ALL, "de_DE");
        }
        tld_strbuf* o = NULL;
        tld_strbuf_alloc(&o, 1024);
        svg_axis xa = {0, 1.5, 0, "x"}, ya = {0, 0.5, 0, "y"};
        svg_grid g;
        svg_panel p;
        double y[3] = {0.25, 0.5, 0.125};
        svg_series s = {"s", "#000000", NULL, NULL, 0, 0.5, y, 3};
        svg_grid_init(&g, 1, 1, 300, 200, 12, 0);
        svg_begin(o, "loc", 300, 200, 12, "t", "d");
        svg_grid_panel(o, &g, 0, 0, NULL, &xa, &ya, &p);
        svg_line(o, &p, &s);
        svg_end(o);
        char err[200];
        CHECK(svg_check(TLD_STR(o), err, sizeof(err)), "%s (locale %s)", err, loc ? loc : "C");
        CHECK(strstr(TLD_STR(o), ">0.5<") != NULL, "tick label 0.5 (locale %s)", loc ? loc : "C");
        setlocale(LC_ALL, "C");
        tld_strbuf_free(o);

        /* 1000-point lines and 1216-bin histograms stay small */
        static double ys[NBIN], hy[1216], he[1217];
        for (int i = 0; i < NBIN; i++) {
                ys[i] = 20.0 + 15.0 * rnd();
        }
        for (int i = 0; i < 1216; i++) {
                hy[i] = 1000.0 * rnd();
        }
        for (int i = 0; i <= 1216; i++) {
                he[i] = 1.0 + i * 0.5;
        }
        tld_strbuf_alloc(&o, 1024);
        svg_axis xb = {1, 1000, 0, NULL}, yb = {0, 40, 0, NULL}, xc = {1, 617, 0, NULL}, yc = {0, 1000, 0, NULL};
        svg_grid_init(&g, 1, 1, 640, 300, 12, 0);
        svg_grid_panel(o, &g, 0, 0, NULL, &xb, &yb, &p);
        size_t before = (size_t)TLD_STRLEN(o);
        svg_series sl = {"l", "#000000", NULL, NULL, 1, 1, ys, NBIN};
        CHECK(svg_line(o, &p, &sl) == OK, "line 1000");
        size_t line_bytes = (size_t)TLD_STRLEN(o) - before;
        CHECK(line_bytes < 13000, "1000-point line is %zu bytes", line_bytes);
        tld_strbuf_clear(o);
        svg_grid_panel(o, &g, 0, 0, NULL, &xc, &yc, &p);
        before = (size_t)TLD_STRLEN(o);
        svg_series sh = {"h", "#000000", NULL, he, 0, 0, hy, 1216};
        CHECK(svg_step(o, &p, &sh) == OK, "step 1216");
        size_t step_bytes = (size_t)TLD_STRLEN(o) - before;
        CHECK(step_bytes < 25000, "1216-bin step histogram is %zu bytes", step_bytes);
        fprintf(stderr, "1000-point line: %zu bytes; 1216-bin step: %zu bytes\n", line_bytes, step_bytes);
        tld_strbuf_free(o);
}

static void test_swatches(void)
{
        tld_strbuf* o = NULL;
        char err[200];
        tld_strbuf_alloc(&o, 1024);
        CHECK(svg_swatch(o, "s1", "#D55E00", "2 3", "MAPQ < 10", "Vermillion dotted line") == OK, "swatch");
        CHECK(svg_check(TLD_STR(o), err, sizeof(err)), "%s", err);
        CHECK(strstr(TLD_STR(o), "stroke-dasharray=\"2 3\"") != NULL && strstr(TLD_STR(o), "stroke=\"#D55E00\"") != NULL, "swatch colour and dash");
        tld_strbuf_clear(o);
        svg_swatch_item it[2] = {{"A", "#009E73", NULL}, {"N", "#2E86C1", "1 2"}};
        CHECK(svg_swatch_row(o, "r1", it, 2, 12, "Bases", "A and N") == OK, "row");
        CHECK(svg_check(TLD_STR(o), err, sizeof(err)), "%s", err);
        CHECK(count_of(TLD_STR(o), "<path") == 2 && count_of(TLD_STR(o), ">A</text>") == 1 && count_of(TLD_STR(o), ">N</text>") == 1, "row items");
        tld_strbuf_clear(o);
        CHECK(svg_swatch_row(o, "r2", NULL, 0, 12, "None", "empty") == OK, "empty row");
        CHECK(svg_check(TLD_STR(o), err, sizeof(err)), "%s", err);
        tld_strbuf_free(o);
}

int main(int argc, char** argv)
{
        test_ticks();
        test_series_paths();
        test_escaping_and_errors();
        test_checker();
        test_grid();
        test_swatches();
        test_tick_max();
        test_specimens(argc > 1 ? argv[1] : NULL);
        test_locale_and_size();
        if (failures) {
                fprintf(stderr, "svg_plot_test: %d check(s) failed\n", failures);
                return EXIT_FAILURE;
        }
        fprintf(stderr, "svg_plot_test: all checks passed\n");
        return EXIT_SUCCESS;
}
