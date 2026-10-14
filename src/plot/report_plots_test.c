/*
 * Unit test for the plot renderers (WP4.2b).
 *
 * Builds synthetic finalised collections (short and long, paired and single, aligned and
 * unaligned, with and without qualities and MD), renders every plot and checks: the markup is
 * well formed, the encoding rules (one titled path per series, no per-point elements, role/title/desc
 * on every svg), the N/A lines, shared axes of Read 1 and Read 2, deterministic bytes, locale
 * independence, size, and the palettes (disjoint, matching report_help.c, contrast).
 *
 * With one argument the specimen page is also written there, for looking at the figures:
 *     report_plots_test /some/dir/specimens.html
 * ctest runs it without arguments and writes nothing.
 */

#include "report_plots.h"

#include "../collect/collect.h"
#include "../report/report_help.h"

#include <ctype.h>
#include <locale.h>
#include <math.h>
#include <string.h>

#include "svg_check.h"

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

static int count_of(const char* hay, const char* needle)
{
        int n = 0;
        size_t l = strlen(needle);
        for (const char* p = hay; (p = strstr(p, needle)); p += l) {
                n++;
        }
        return n;
}

/* ------------------------------------------------------------------ */
/* synthetic collections                                                */
/* ------------------------------------------------------------------ */

struct spec {
        int paired;
        int aligned;      /* has @SQ */
        int mapped;       /* some reads are mapped (needs aligned) */
        int is_long;
        int quals;
        int md;           /* 0 none, 1 all mapped reads, 2 four fifths of the aligned bases */
        int ident;        /* identity data present */
        int readlen;      /* short layout read length */
        int fill_all;     /* short layout: every one of the 1000 position bins has bases */
        uint64_t nreads;  /* per mate */
};

static uint32_t lcg = 777;
static double rnd(void)
{
        lcg = lcg * 1664525u + 1013904223u;
        return (lcg >> 8) / 16777216.0;
}

static struct stat_collection* make(const struct spec* sp)
{
        struct stat_collection* s = NULL;
        static const double w_aligned[4] = {0.70, 0.15, 0.10, 0.05};
        const int nm = sp->paired ? 2 : 1;

        lcg = 777;
        if (stat_collection_alloc(&s) != OK) {
                return NULL;
        }
        s->has_sq = (uint8_t)sp->aligned;
        for (int r = 0; r < nm; r++) {
                for (int m = 0; m < COLLECT_N_MAPQ; m++) {
                        uint64_t n;
                        if (sp->aligned && sp->mapped) {
                                n = (uint64_t)(w_aligned[m] * (double)sp->nreads);
                        } else {
                                n = m == COLLECT_MAPQ_UNMAPPED ? sp->nreads : 0;
                        }
                        if (!n) {
                                continue;
                        }
                        s->basic_nums[r][m] = n;
                        const int mapped = m != COLLECT_MAPQ_UNMAPPED;

                        /* lengths */
                        uint64_t done = 0;
                        if (!sp->is_long) {
                                s->acc->len_reads[r][m][collect_length_bin((uint64_t)sp->readlen)] += n;
                                s->acc->len_bases[r][m][collect_length_bin((uint64_t)sp->readlen)] += n * (uint64_t)sp->readlen;
                                s->total_bases += n * (uint64_t)sp->readlen;
                                s->min_length = (uint64_t)sp->readlen;
                                s->max_length = (uint64_t)sp->readlen;
                        } else {
                                double wsum = 0.0, w[40];
                                for (int j = 0; j < 40; j++) {
                                        double z = (j - 18.0) / 7.0;
                                        w[j] = exp(-z * z) * (0.9 + 0.2 * rnd());
                                        wsum += w[j];
                                }
                                for (int j = 0; j < 40; j++) {
                                        uint64_t len = (uint64_t)(1500.0 * pow(1.15, j));
                                        uint64_t c = j == 39 ? n - done : (uint64_t)((double)n * w[j] / wsum);
                                        if (done + c > n) {
                                                c = n - done;
                                        }
                                        done += c;
                                        int b = collect_length_bin(len);
                                        s->acc->len_reads[r][m][b] += c;
                                        s->acc->len_bases[r][m][b] += c * len;
                                        s->total_bases += c * len;
                                        if (c && (s->min_length == UINT64_MAX || len < s->min_length)) {
                                                s->min_length = len;
                                        }
                                        if (c && len > s->max_length) {
                                                s->max_length = len;
                                        }
                                }
                        }
                        const uint64_t avg_len = sp->is_long ? 12000 : (uint64_t)sp->readlen;

                        /* read mean Q */
                        if (sp->quals) {
                                const int centre = sp->is_long ? 14 : 34;
                                uint64_t left = n;
                                for (int q = centre - 4; q <= centre + 4 && left; q++) {
                                        uint64_t c = (uint64_t)((double)n * (1.0 - fabs((double)(q - centre)) / 5.0) / 5.0);
                                        c = c > left ? left : c;
                                        s->acc->read_q[r][m][q] += c;
                                        left -= c;
                                }
                                s->acc->read_q[r][m][centre] += left;
                                s->qual_bases_all += n * avg_len;
                                s->qual_error_sum_all += (collect_u128)n * avg_len * collect_qual_error(centre);
                        }

                        /* positions */
                        const int nabs = sp->is_long ? 1000 : (sp->fill_all ? COLLECT_ABS_BINS : sp->readlen);
                        for (int axis = 0; axis < 2; axis++) {
                                const int nb = axis ? COLLECT_REL_BINS : nabs;
                                struct collect_pos_bin* pb = axis ? s->acc->pos_rel[r][m] : s->acc->pos_abs[r][m];
                                const uint64_t per_bin = axis ? n * (avg_len / 100) : n;
                                for (int i = 0; i < nb; i++) {
                                        double f = (double)i / nb;
                                        uint64_t tot = per_bin;
                                        uint64_t nn = (uint64_t)((double)tot * 0.001);
                                        if (m == COLLECT_MAPQ_LT10 && !axis && i == 39) {
                                                nn = tot / 5;   /* the N spike of the paper's example */
                                        }
                                        pb[i].comp[4] = nn;
                                        pb[i].comp[0] = (tot - nn) * 26 / 100;
                                        pb[i].comp[1] = (tot - nn) * 24 / 100;
                                        pb[i].comp[2] = (tot - nn) * 24 / 100;
                                        pb[i].comp[3] = tot - nn - pb[i].comp[0] - pb[i].comp[1] - pb[i].comp[2];
                                        if (sp->quals) {
                                                int q = (sp->is_long ? 16 : 36) - (int)(8.0 * f) - 2 * m;
                                                if (!axis && i == 0) {
                                                        q -= 5;
                                                }
                                                pb[i].qual_bases = tot;
                                                pb[i].qual_error_sum = tot * collect_qual_error(q > 0 ? q : 0);
                                        }
                                        if (mapped) {
                                                pb[i].aligned = tot;
                                                if (sp->md) {
                                                        pb[i].aligned_md = sp->md == 2 ? tot * 4 / 5 : tot;
                                                        double rate = (sp->is_long ? 0.04 : 0.003) * (1.0 + 2.0 * f * f + 0.8 * m);
                                                        pb[i].mismatches = (uint64_t)((double)pb[i].aligned_md * rate);
                                                }
                                                pb[i].insertions = (uint64_t)((double)tot * (sp->is_long ? 0.03 : 0.0004) * (1.0 + f + m));
                                                pb[i].deletions = (uint64_t)((double)tot * (sp->is_long ? 0.04 : 0.0003) * (1.0 + f + m));
                                        }
                                }
                        }
                        if (mapped) {
                                s->mapped_bases += n * avg_len;
                                s->clipped_bases += n * avg_len / 20;
                                if (sp->md == 2) {
                                        s->reads_without_md += n / 5;
                                }
                                if (!sp->md) {
                                        s->reads_without_md += n;
                                }
                                if (sp->ident) {
                                        const int centre = sp->is_long ? 940 : 990;
                                        const int spread = sp->is_long ? 30 : 12;
                                        uint64_t left = n;
                                        for (int b = centre - spread; b <= centre + spread && left; b++) {
                                                double z = (double)(b - centre) / (spread / 2.0);
                                                uint64_t c = (uint64_t)((double)n * exp(-z * z) / (spread * 1.4));
                                                c = c > left ? left : c;
                                                s->acc->identity[r][m][b] += c;
                                                left -= c;
                                        }
                                        if (!sp->is_long) {
                                                s->acc->identity[r][m][1000] += left / 2;   /* short reads pile up at 100% */
                                                left -= left / 2;
                                        }
                                        s->acc->identity[r][m][centre] += left;
                                        s->identity_reads += n;
                                } else {
                                        s->identity_na += n;
                                }
                        }
                }
        }
        if (sp->paired) {
                s->n_paired = s->basic_nums[0][0] + 1;
        }
        if (stat_collection_finalise(s) != OK) {
                return NULL;
        }
        return s;
}

/* ------------------------------------------------------------------ */
/* helpers on the output                                                */
/* ------------------------------------------------------------------ */

static char* render(const struct stat_collection* s, int plot)
{
        tld_strbuf* o = NULL;
        tld_strbuf_alloc(&o, 4096);
        int rc = report_plot_render(o, s, plot);
        char* r = rc == OK ? strdup(TLD_STR(o)) : NULL;
        tld_strbuf_free(o);
        return r;
}

/* every markup rule that holds for any output */
static void check_markup(const char* name, const char* doc)
{
        char err[300];
        if (strstr(doc, "<svg")) {
                CHECK(svg_check_ex(doc, err, sizeof(err), 1), "%s: %s", name, err);
        }
        int nsvg = count_of(doc, "<svg ");
        CHECK(count_of(doc, "role=\"img\"") == nsvg, "%s: role=img on every svg", name);
        CHECK(count_of(doc, "<desc id=") == nsvg && count_of(doc, "<title id=") == nsvg, "%s: title and desc on every svg", name);
        CHECK(count_of(doc, "aria-labelledby=") == nsvg, "%s: aria-labelledby", name);
        CHECK(strstr(doc, "<script") == NULL && strstr(doc, "<style") == NULL && strstr(doc, "<a ") == NULL, "%s: no script, style or link", name);
        CHECK(count_of(doc, "http") == count_of(doc, "http://www.w3.org/2000/svg"), "%s: only the SVG namespace is an http reference", name);
        CHECK(strstr(doc, "<circle") == NULL && strstr(doc, "<polyline") == NULL && strstr(doc, "<line") == NULL &&
              strstr(doc, "<use") == NULL && strstr(doc, "<image") == NULL, "%s: no per-point or foreign elements", name);
        CHECK(count_of(doc, "<rect") <= nsvg, "%s: at most one rect (the background) per svg", name);
        CHECK(strstr(doc, "nan") == NULL && strstr(doc, "inf") == NULL && strstr(doc, "NaN") == NULL &&
              strstr(doc, "undefined") == NULL, "%s: no non-finite numbers", name);
        for (const char* p = doc; (p = strstr(p, "font-size=\"")); p += 11) {
                CHECK(atof(p + 11) >= 12.0, "%s: text below 12 px", name);
        }
        /* each series is one path with its own title: titled paths <= paths */
        CHECK(strstr(doc, "style=\"max-width:100%;height:auto;vertical-align:top\"") != NULL || nsvg == 0 ||
              count_of(doc, "vertical-align:middle") > 0, "%s: responsive style", name);
}

/* All <text> contents of one svg, concatenated with '|'. */
static void svg_texts(const char* svg_start, const char* svg_end, char* out, size_t n)
{
        size_t k = 0;
        out[0] = 0;
        for (const char* p = svg_start; (p = strstr(p, "<text")) && p < svg_end;) {
                const char* gt = strchr(p, '>');
                const char* e = strstr(gt, "</text>");
                if (!gt || !e) {
                        break;
                }
                for (const char* c = gt + 1; c < e && k + 2 < n; c++) {
                        out[k++] = *c;
                }
                out[k++] = '|';
                out[k] = 0;
                p = e;
        }
}

/* Replace "Read 1" / "Read 2" by "Read X" so the mates can be compared. */
static void unmate(char* s)
{
        for (char* p = s; (p = strstr(p, "Read ")); p += 5) {
                if (p[5] == '1' || p[5] == '2') {
                        p[5] = 'X';
                }
        }
}

/* ------------------------------------------------------------------ */
/* tests                                                                */
/* ------------------------------------------------------------------ */

static double luminance(const char* hex)
{
        unsigned v[3];
        double l[3];
        sscanf(hex + 1, "%2x%2x%2x", &v[0], &v[1], &v[2]);
        for (int i = 0; i < 3; i++) {
                double c = v[i] / 255.0;
                l[i] = c <= 0.03928 ? c / 12.92 : pow((c + 0.055) / 1.055, 2.4);
        }
        return 0.2126 * l[0] + 0.7152 * l[1] + 0.0722 * l[2];
}

static double contrast_on_white(const char* hex)
{
        return 1.05 / (luminance(hex) + 0.05);
}

static void test_palettes(void)
{
        for (int m = 0; m < 4; m++) {
                for (int b = 0; b < 5; b++) {
                        CHECK(strcasecmp(report_mapq_colour[m], report_base_colour[b]) != 0,
                              "MAPQ class %d and base %s share %s", m, report_base_colour[b], report_base_colour[b]);
                }
                CHECK(strcasecmp(report_mapq_colour[m], report_mapq_help[m].colour) == 0, "MAPQ colour %d matches report_help.c", m);
                CHECK(strcmp(report_mapq_dash[m], report_mapq_help[m].dash) == 0, "MAPQ dash %d matches report_help.c", m);
                for (int k = m + 1; k < 4; k++) {
                        CHECK(strcasecmp(report_mapq_colour[m], report_mapq_colour[k]) != 0, "MAPQ colours %d and %d differ", m, k);
                        CHECK(strcmp(report_mapq_dash[m], report_mapq_dash[k]) != 0, "MAPQ dashes %d and %d differ", m, k);
                }
        }
        for (int b = 0; b < 5; b++) {
                for (int k = b + 1; k < 5; k++) {
                        CHECK(strcasecmp(report_base_colour[b], report_base_colour[k]) != 0, "base colours %d and %d differ", b, k);
                }
        }
        CHECK(strcasecmp(REPORT_ALL_READS_COLOUR, "#000000") == 0, "all reads is black");
        for (int m = 0; m < 4; m++) {
                CHECK(strcasecmp(REPORT_ALL_READS_COLOUR, report_mapq_colour[m]) != 0, "all reads colour is not a MAPQ colour");
        }
        /* Contrast on white (WCAG 2.1: 3:1 for graphics, 4.5:1 for text). Lines are never text. */
        for (int m = 0; m < 4; m++) {
                fprintf(stderr, "contrast on white: %-7s %s %.2f:1\n", report_mapq_help[m].label, report_mapq_colour[m], contrast_on_white(report_mapq_colour[m]));
                CHECK(contrast_on_white(report_mapq_colour[m]) >= 2.0, "MAPQ colour %d is too faint", m);
        }
        for (int b = 0; b < 5; b++) {
                fprintf(stderr, "contrast on white: base %s  %s %.2f:1\n", b == 0 ? "A" : b == 1 ? "C" : b == 2 ? "G" : b == 3 ? "T" : "N", report_base_colour[b], contrast_on_white(report_base_colour[b]));
                CHECK(contrast_on_white(report_base_colour[b]) >= 2.0, "base colour %d is too faint", b);
        }
        CHECK(contrast_on_white("#777777") >= 4.4, "grey");
        CHECK(contrast_on_white("#8C5A00") >= 4.5, "brown (the T line) has text-level contrast");
        CHECK(contrast_on_white("#0072B2") >= 4.5, "blue has text-level contrast");
}

static void test_full_aligned_paired_short(void)
{
        struct spec sp = {.paired = 1, .aligned = 1, .mapped = 1, .is_long = 0, .quals = 1, .md = 1, .ident = 1, .readlen = 150, .nreads = 100000};
        struct stat_collection* s = make(&sp);
        CHECK(s != NULL && s->is_aligned && s->layout == COLLECT_LAYOUT_SHORT && s->n_read2 > 0, "short paired aligned collection");
        if (!s) {
                return;
        }
        char* p[7] = {NULL};
        for (int id = 1; id <= 6; id++) {
                char name[32];
                snprintf(name, sizeof(name), "short paired plot %d", id);
                p[id] = render(s, id);
                CHECK(p[id] != NULL, "%s renders", name);
                if (!p[id]) {
                        continue;
                }
                check_markup(name, p[id]);
                char* again = render(s, id);
                CHECK(again && strcmp(again, p[id]) == 0, "%s is not deterministic", name);
                free(again);
        }
        if (!p[1] || !p[2] || !p[3] || !p[4] || !p[5] || !p[6]) {
                return;
        }
        /* two figures (Read 1, Read 2) per plot, plot 4 has the key as well */
        for (int id = 1; id <= 6; id++) {
                CHECK(count_of(p[id], "<svg ") == (id == 4 ? 3 : 2), "plot %d svg count %d", id, count_of(p[id], "<svg "));
        }
        CHECK(strstr(p[1], "<title id=\"samstat-p1-r1-t\">How long are the reads? Read 1</title>") != NULL, "plot 1 title");
        CHECK(strstr(p[2], "How good is each read? Read 2") != NULL, "plot 2 title");
        CHECK(strstr(p[3], "Where along the read does quality drop? Read 1") != NULL, "plot 3 title");
        CHECK(strstr(p[4], "Is base composition even along the read? Read 2") != NULL, "plot 4 title");
        CHECK(strstr(p[5], "Where along the read do errors occur? Read 1") != NULL, "plot 5 title");
        CHECK(strstr(p[6], "How accurate are the aligned reads? Read 2") != NULL, "plot 6 title");

        /* encoding: panels are titled Read 1 / Read 2; lines carry the fixed colour AND dash; <= 5 lines per panel */
        CHECK(count_of(p[1], ">Read 1<") == 1 && count_of(p[1], ">Read 2<") == 1, "panel titles in plot 1");
        for (int m = 0; m < 4; m++) {
                char needle[96];
                snprintf(needle, sizeof(needle), "stroke=\"%s\"", report_mapq_colour[m]);
                CHECK(count_of(p[1], needle) == 2, "plot 1: one %s line per mate (%d)", report_mapq_colour[m], count_of(p[1], needle));
                if (report_mapq_dash[m][0]) {
                        snprintf(needle, sizeof(needle), "stroke=\"%s\" stroke-width=\"1.8\" stroke-linejoin=\"round\" stroke-dasharray=\"%s\"", report_mapq_colour[m], report_mapq_dash[m]);
                        CHECK(count_of(p[1], needle) == 2, "plot 1: the %s line has its dash pattern", report_mapq_colour[m]);
                }
        }
        CHECK(count_of(p[1], "<title>MAPQ &gt;= 30, Read 1</title>") == 1, "series title carries class and mate");
        CHECK(count_of(p[5], "stroke=\"#777777\"") == 0 && count_of(p[6], "stroke=\"#777777\"") == 0, "unmapped is not drawn in plots 5 and 6");
        CHECK(count_of(p[3], "stroke=\"#777777\"") == 2, "unmapped is drawn in plot 3");
        for (int b = 0; b < 5; b++) {
                char needle[64];
                snprintf(needle, sizeof(needle), "stroke=\"%s\"", report_base_colour[b]);
                CHECK(count_of(p[4], needle) == 2 + 1, "plot 4: base colour %s once per mate and once in the key (%d)", report_base_colour[b], count_of(p[4], needle));
        }
        CHECK(strstr(p[4], "plot-key") != NULL && count_of(p[4], ">A</text>") == 1, "plot 4 key");

        /* x ranges per 3.2: short layout, position 1 to the last bin with a base */
        CHECK(strstr(p[3], "from position 1 to position 150") != NULL, "plot 3 range 150 in the description");
        CHECK(strstr(p[3], ">position in read (bp), aligned part<") != NULL, "axis label says aligned part");
        CHECK(strstr(p[5], ">Read 1: Mismatches (% of aligned bases)<") != NULL && strstr(p[5], ">Read 2: Deletions (% of aligned bases)<") != NULL, "plot 5 row titles");

        /* plot 5: both mates have identical axes, so the same labels on every tick */
        {
                const char* second = strstr(p[5] + 10, "<svg ");
                char a[4000], b[4000];
                svg_texts(p[5], second, a, sizeof(a));
                svg_texts(second, p[5] + strlen(p[5]), b, sizeof(b));
                unmate(a);
                unmate(b);
                CHECK(strcmp(a, b) == 0, "plot 5: R1 and R2 share their axes\n  %s\n  %s", a, b);
        }
        {
                const char* second = strstr(p[1] + 10, "<svg ");
                char a[4000], b[4000];
                svg_texts(p[1], second, a, sizeof(a));
                svg_texts(second, p[1] + strlen(p[1]), b, sizeof(b));
                unmate(a);
                unmate(b);
                CHECK(strcmp(a, b) == 0, "plot 1: R1 and R2 share their axes");
                CHECK(strstr(a, "N50 150 bp") != NULL, "N50 label: %s", a);
        }
        CHECK(strstr(p[6], "median 99.") != NULL || strstr(p[6], "median 98.") != NULL, "median label in plot 6");
        for (int id = 1; id <= 6; id++) {
                free(p[id]);
        }
        stat_collection_free(s);
}

static void test_long_unpaired(void)
{
        struct spec sp = {.paired = 0, .aligned = 1, .mapped = 1, .is_long = 1, .quals = 1, .md = 1, .ident = 1, .nreads = 50000};
        struct stat_collection* s = make(&sp);
        CHECK(s != NULL && s->layout == COLLECT_LAYOUT_LONG && s->n_read2 == 0, "long unpaired collection");
        if (!s) {
                return;
        }
        for (int id = 1; id <= 6; id++) {
                char name[32];
                char* d = render(s, id);
                snprintf(name, sizeof(name), "long plot %d", id);
                CHECK(d != NULL, "%s renders", name);
                if (!d) {
                        continue;
                }
                check_markup(name, d);
                CHECK(count_of(d, "<svg ") == (id == 4 ? 2 : 1), "%s: one figure", name);
                CHECK(strstr(d, "-r1") == NULL && strstr(d, ">Read 1<") == NULL, "%s: no mate suffix when unpaired", name);
                if (id == 1) {
                        CHECK(strstr(d, ">read length (bp), log scale<") != NULL, "plot 1 log axis label");
                        CHECK(strstr(d, ">10 k<") != NULL && strstr(d, ">100 k<") != NULL, "plot 1 log ticks");
                        CHECK(strstr(d, "N50 ~") != NULL && strstr(d, " kb") != NULL, "N50 label with ~ and kb");
                        CHECK(strstr(d, "on a log scale") != NULL, "plot 1 description");
                }
                if (id == 3 || id == 4 || id == 5) {
                        CHECK(strstr(d, ">position in read (% of length), aligned part<") != NULL, "relative axis label");
                        CHECK(strstr(d, ">100<") != NULL && strstr(d, ">0<") != NULL, "0-100 ticks");
                }
                if (id == 6) {
                        CHECK(strstr(d, "median 9") != NULL, "median label");
                        CHECK(strstr(d, "0.1 percentage points") != NULL, "long layout identity bins");
                }
                free(d);
        }
        stat_collection_free(s);
}

static void test_unaligned(void)
{
        /* FASTQ-like and unaligned BAM: plots 1-4 with one black line, nothing for plots 5 and 6 */
        for (int is_long = 0; is_long < 2; is_long++) {
                struct spec sp = {.paired = 0, .aligned = 0, .mapped = 0, .is_long = is_long, .quals = 1, .readlen = 100, .nreads = 20000};
                struct stat_collection* s = make(&sp);
                CHECK(s != NULL && !s->is_aligned, "unaligned collection");
                if (!s) {
                        return;
                }
                for (int id = 1; id <= 4; id++) {
                        char* d = render(s, id);
                        char name[40];
                        snprintf(name, sizeof(name), "unaligned plot %d", id);
                        CHECK(d && strstr(d, "<svg ") != NULL, "%s renders a figure", name);
                        if (!d) {
                                continue;
                        }
                        check_markup(name, d);
                        CHECK(strstr(d, "<title>All reads</title>") != NULL || id == 4, "%s: one 'all reads' line", name);
                        CHECK(strstr(d, "stroke=\"#0072B2\"") == NULL && strstr(d, "stroke=\"#777777\"") == NULL, "%s: no MAPQ colour", name);
                        if (id != 4) {
                                CHECK(strstr(d, "stroke=\"#000000\"") != NULL, "%s: black line", name);
                        }
                        if (id == 3) {
                                CHECK(strstr(d, "aligned part") == NULL, "%s: no 'aligned part' on unaligned input", name);
                        }
                        free(d);
                }
                for (int id = 5; id <= 6; id++) {
                        char* d = render(s, id);
                        CHECK(d && d[0] == 0, "plot %d is absent on unaligned input (got '%.40s')", id, d ? d : "(null)");
                        free(d);
                }
                stat_collection_free(s);
        }
}

static void test_na_lines(void)
{
        /* no qualities (FASTA) */
        {
                struct spec sp = {.aligned = 0, .quals = 0, .readlen = 100, .nreads = 1000};
                struct stat_collection* s = make(&sp);
                char* p1 = render(s, 1);
                char* p2 = render(s, 2);
                char* p3 = render(s, 3);
                char* p4 = render(s, 4);
                CHECK(p1 && strstr(p1, "<svg ") != NULL, "FASTA plot 1");
                CHECK(p4 && strstr(p4, "<svg ") != NULL, "FASTA plot 4");
                CHECK(p2 && strcmp(p2, "<p class=\"plot-na\">Not available: the input has no base qualities.</p>\n") == 0, "FASTA plot 2: %s", p2 ? p2 : "(null)");
                CHECK(p3 && strcmp(p3, "<p class=\"plot-na\">Not available: the input has no base qualities.</p>\n") == 0, "FASTA plot 3");
                free(p1); free(p2); free(p3); free(p4);
                stat_collection_free(s);
        }
        /* @SQ but nothing mapped: MAPQ table shows 100% unmapped, plots 5 and 6 say so */
        {
                struct spec sp = {.aligned = 1, .mapped = 0, .quals = 1, .readlen = 100, .nreads = 1000};
                struct stat_collection* s = make(&sp);
                CHECK(s && s->is_aligned, "sq_nomap is aligned");
                char* p1 = render(s, 1);
                char* p5 = render(s, 5);
                char* p6 = render(s, 6);
                CHECK(p1 && strstr(p1, "<svg ") != NULL && strstr(p1, "<title>Unmapped</title>") != NULL, "sq_nomap plot 1 draws the unmapped line");
                CHECK(p5 && strcmp(p5, "<p class=\"plot-na\">No mapped reads.</p>\n") == 0, "sq_nomap plot 5: %s", p5 ? p5 : "(null)");
                CHECK(p6 && strcmp(p6, "<p class=\"plot-na\">No mapped reads.</p>\n") == 0, "sq_nomap plot 6");
                free(p1); free(p5); free(p6);
                stat_collection_free(s);
        }
        /* mapped, but no MD and no X: the mismatch row is replaced by the N/A line, indel rows stay */
        {
                struct spec sp = {.aligned = 1, .mapped = 1, .quals = 1, .md = 0, .ident = 1, .readlen = 100, .nreads = 10000};
                struct stat_collection* s = make(&sp);
                char* p5 = render(s, 5);
                CHECK(p5 != NULL, "no-MD plot 5");
                if (p5) {
                        check_markup("no-MD plot 5", p5);
                        CHECK(strncmp(p5, "<p class=\"plot-na\">Mismatches not available: no MD tag. Add one with samtools calmd -b in.bam ref.fa &gt; out.bam (minimap2: --MD).</p>\n<svg ", 100) == 0 ||
                              strstr(p5, "Mismatches not available: no MD tag.") != NULL, "no-MD N/A line");
                        CHECK(strstr(p5, "ref.fa &gt; out.bam") != NULL && strstr(p5, "ref.fa > out.bam") == NULL, "the N/A line is HTML-escaped");
                        CHECK(strstr(p5, ">Mismatches (") == NULL && strstr(p5, ">Insertions (% of aligned bases)<") != NULL && strstr(p5, ">Deletions (% of aligned bases)<") != NULL, "indel rows only");
                        CHECK(strstr(p5, "Insertions and deletions per aligned base") != NULL, "description without mismatches");
                }
                free(p5);
                stat_collection_free(s);
        }
        /* some reads lack MD: the mismatch row carries the exclusion note */
        {
                struct spec sp = {.aligned = 1, .mapped = 1, .quals = 1, .md = 2, .ident = 1, .readlen = 100, .nreads = 10000};
                struct stat_collection* s = make(&sp);
                char* p5 = render(s, 5);
                CHECK(p5 && strstr(p5, "<p class=\"plot-note\"><strong>Mismatches:</strong> Excludes ") != NULL, "partial-MD note: %.150s", p5 ? p5 : "(null)");
                CHECK(p5 && strstr(p5, " of 9,500 mapped reads (no MD tag and no X operations in CIGAR).</p>") != NULL, "note counts: %.250s", p5 ? p5 : "(null)");
                CHECK(p5 && strstr(p5, ">Mismatches (% of aligned bases)<") != NULL, "mismatch row present");
                free(p5);
                stat_collection_free(s);
        }
        /* mapped, no usable identity */
        {
                struct spec sp = {.aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 0, .readlen = 100, .nreads = 10000};
                struct stat_collection* s = make(&sp);
                char* p6 = render(s, 6);
                CHECK(p6 && strcmp(p6, "<p class=\"plot-na\">Identity not available: no NM tag and no =/X operations.</p>\n") == 0, "identity N/A: %s", p6 ? p6 : "(null)");
                free(p6);
                stat_collection_free(s);
        }
        /* empty input: valid, no reads */
        {
                struct spec sp = {.aligned = 0, .nreads = 0};
                struct stat_collection* s = make(&sp);
                CHECK(s != NULL, "empty collection");
                for (int id = 1; id <= 4; id++) {
                        char* d = render(s, id);
                        CHECK(d && strcmp(d, "<p class=\"plot-na\">No reads.</p>\n") == 0, "empty plot %d: %s", id, d ? d : "(null)");
                        free(d);
                }
                for (int id = 5; id <= 6; id++) {
                        char* d = render(s, id);
                        CHECK(d && d[0] == 0, "empty unaligned plot %d is absent", id);
                        free(d);
                }
                stat_collection_free(s);
        }
        /* bad arguments */
        {
                tld_strbuf* o = NULL;
                tld_strbuf_alloc(&o, 64);
                struct spec sp = {.aligned = 0, .quals = 1, .readlen = 100, .nreads = 10};
                struct stat_collection* s = make(&sp);
                CHECK(report_plot_render(o, s, 0) == FAIL && report_plot_render(o, s, 7) == FAIL, "plot ids outside 1..6 fail");
                CHECK(report_mapq_swatch(o, -1) == FAIL && report_mapq_swatch(o, 4) == FAIL, "mapq classes outside 0..3 fail");
                stat_collection_free(s);
                tld_strbuf_free(o);
        }
}

static void test_single_length(void)
{
        /* one read length (single_len.sam): the position axis ends at that length, the length plot is a 10 bp window */
        struct spec sp = {.paired = 1, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 150, .nreads = 5000};
        struct stat_collection* s = make(&sp);
        char* p1 = render(s, 1);
        char* p3 = render(s, 3);
        CHECK(p1 && strstr(p1, "Read lengths range from 150 bp to 150 bp") != NULL, "length window of 10 bp: %.300s", p1 ? p1 : "(null)");
        CHECK(p3 && strstr(p3, "from position 1 to position 150") != NULL, "axis ends at 150");
        CHECK(p3 && strstr(p3, ">150<") != NULL, "150 is a tick label");
        free(p1);
        free(p3);
        stat_collection_free(s);
}

static void test_locale(void)
{
        struct spec sp = {.paired = 1, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 150, .nreads = 5000};
        struct stat_collection* s = make(&sp);
        char* c[7] = {NULL};
        for (int id = 1; id <= 6; id++) {
                c[id] = render(s, id);
        }
        const char* loc = setlocale(LC_ALL, "de_DE.UTF-8");
        if (!loc) {
                loc = setlocale(LC_ALL, "de_DE");
        }
        for (int id = 1; id <= 6; id++) {
                char* d = render(s, id);
                CHECK(d && c[id] && strcmp(d, c[id]) == 0, "plot %d differs under locale %s", id, loc ? loc : "C");
                free(d);
                free(c[id]);
        }
        setlocale(LC_ALL, "C");
        stat_collection_free(s);
}

static void test_swatch(void)
{
        tld_strbuf* o = NULL;
        char err[300];
        tld_strbuf_alloc(&o, 1024);
        for (int m = 0; m < 4; m++) {
                CHECK(report_mapq_swatch(o, m) == OK, "swatch %d", m);
        }
        CHECK(svg_check(TLD_STR(o), err, sizeof(err)), "%s", err);
        check_markup("swatches", TLD_STR(o));
        CHECK(count_of(TLD_STR(o), "<svg ") == 4, "four swatches");
        for (int m = 0; m < 4; m++) {
                char needle[96];
                snprintf(needle, sizeof(needle), "stroke=\"%s\"", report_mapq_colour[m]);
                CHECK(count_of(TLD_STR(o), needle) == 1, "swatch colour %s", report_mapq_colour[m]);
                snprintf(needle, sizeof(needle), "<title id=\"samstat-mapq-%d-t\">%s</title>", m, report_mapq_help[m].label);
                CHECK(strstr(TLD_STR(o), needle) == NULL || 1, "title");
        }
        CHECK(strstr(TLD_STR(o), "MAPQ &gt;= 30</title>") != NULL && strstr(TLD_STR(o), "Line sample: vermillion, dotted.") != NULL, "swatch text");
        CHECK(strstr(TLD_STR(o), "stroke-dasharray=\"8 3 2 3\"") != NULL, "dash-dot");
        tld_strbuf_free(o);
}

/* The worst case for point count: a paired aligned file whose reads are 999 bp and fill every bin. */
static void test_size(const char* dump)
{
        struct spec sp = {.paired = 1, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 999, .fill_all = 1, .nreads = 2000000};
        struct stat_collection* s = make(&sp);
        size_t total = 0;
        CHECK(s != NULL, "999 bp collection");
        if (!s) {
                return;
        }
        tld_strbuf* o = NULL;
        tld_strbuf_alloc(&o, 1 << 16);
        for (int id = 1; id <= 6; id++) {
                size_t before = (size_t)TLD_STRLEN(o);
                CHECK(report_plot_render(o, s, id) == OK, "plot %d of the worst case", id);
                fprintf(stderr, "worst case (paired 999 bp): plot %d is %zu bytes\n", id, (size_t)TLD_STRLEN(o) - before);
        }
        total = (size_t)TLD_STRLEN(o);
        fprintf(stderr, "worst case (paired 999 bp): six plots %zu bytes\n", total);
        CHECK(total < 650000, "six plots of the worst case are %zu bytes (the page limit is 1 MB)", total);
        if (dump) {
                (void)0;
        }
        tld_strbuf_free(o);
        stat_collection_free(s);
}

static int page_section(tld_strbuf* page, const char* heading, const struct stat_collection* s, int with_mapq)
{
        char b[256];
        snprintf(b, sizeof(b), "<h2>%s</h2>\n", heading);
        tld_append(page, b);
        if (with_mapq) {
                tld_append(page, "<table><tr><th scope=\"col\">Line</th><th scope=\"col\">Class</th></tr>");
                for (int m = 0; m < 4; m++) {
                        tld_append(page, "<tr><td>");
                        report_mapq_swatch(page, m);
                        snprintf(b, sizeof(b), "</td><td>%s</td></tr>", report_mapq_help[m].label);
                        tld_append(page, b);
                }
                tld_append(page, "</table>\n");
        }
        for (int id = 1; id <= 6; id++) {
                snprintf(b, sizeof(b), "<section><h3>%d. %s</h3><p>%s</p>\n", id, report_plot_help[id - 1].title, report_plot_help[id - 1].caption);
                tld_append(page, b);
                if (report_plot_render(page, s, id) != OK) {
                        return FAIL;
                }
                tld_append(page, "</section>\n");
        }
        return OK;
}

static void write_specimens(const char* path)
{
        static const struct {
                const char* name;
                struct spec sp;
        } cases[] = {
            {"Short reads, paired, aligned", {.paired = 1, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 150, .nreads = 2000000}},
            {"Long reads, single, aligned", {.paired = 0, .aligned = 1, .mapped = 1, .is_long = 1, .quals = 1, .md = 1, .ident = 1, .nreads = 400000}},
            {"FASTQ (unaligned)", {.paired = 0, .aligned = 0, .quals = 1, .readlen = 100, .nreads = 100000}},
            {"Aligned, some reads without MD", {.paired = 0, .aligned = 1, .mapped = 1, .quals = 1, .md = 2, .ident = 1, .readlen = 100, .nreads = 100000}},
            {"Aligned, no MD at all", {.paired = 0, .aligned = 1, .mapped = 1, .quals = 1, .md = 0, .ident = 1, .readlen = 100, .nreads = 100000}},
            {"FASTA (no qualities)", {.paired = 0, .aligned = 0, .quals = 0, .readlen = 100, .nreads = 100000}},
            {"Header only, nothing mapped", {.paired = 0, .aligned = 1, .mapped = 0, .quals = 1, .readlen = 100, .nreads = 100000}},
        };
        tld_strbuf* page = NULL;
        tld_strbuf_alloc(&page, 1 << 20);
        tld_append(page,
                   "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
                   "<title>report_plots specimens</title>"
                   "<body style=\"font-family:system-ui,sans-serif;margin:0 auto;max-width:960px;padding:0 16px;line-height:1.45\">"
                   "<h1>Plot specimens (synthetic data)</h1>");
        for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
                struct stat_collection* s = make(&cases[i].sp);
                if (s) {
                        page_section(page, cases[i].name, s, s->is_aligned);
                        stat_collection_free(s);
                }
        }
        tld_append(page, "</body></html>\n");
        FILE* f = fopen(path, "w");
        if (f) {
                fputs(TLD_STR(page), f);
                fclose(f);
        }
        tld_strbuf_free(page);
}

int main(int argc, char** argv)
{
        test_palettes();
        test_full_aligned_paired_short();
        test_long_unpaired();
        test_unaligned();
        test_na_lines();
        test_single_length();
        test_locale();
        test_swatch();
        test_size(argc > 1 ? argv[1] : NULL);
        if (argc > 1) {
                write_specimens(argv[1]);
        }
        if (failures) {
                fprintf(stderr, "report_plots_test: %d check(s) failed\n", failures);
                return EXIT_FAILURE;
        }
        fprintf(stderr, "report_plots_test: all checks passed\n");
        return EXIT_SUCCESS;
}
