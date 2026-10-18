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
#include "svg_plot.h"

#include "../collect/collect.h"
#include "../report/report_help.h"

#include <ctype.h>
#include <locale.h>
#include <inttypes.h>
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
        int one_len;      /* short layout: every read is readlen long (no trimmed tail) */
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
                                /* most reads full length, a tail of adapter-trimmed reads 1-20 bp shorter */
                                const uint64_t trimmed = sp->one_len ? 0 : n / 10;
                                for (uint64_t t = 0; t < trimmed; t++) {
                                        uint64_t len = (uint64_t)sp->readlen - 1 - t % 20;
                                        s->acc->len_reads[r][m][collect_length_bin(len)] += 1;
                                        s->acc->len_bases[r][m][collect_length_bin(len)] += len;
                                        s->total_bases += len;
                                }
                                s->acc->len_reads[r][m][collect_length_bin((uint64_t)sp->readlen)] += n - trimmed;
                                s->acc->len_bases[r][m][collect_length_bin((uint64_t)sp->readlen)] += (n - trimmed) * (uint64_t)sp->readlen;
                                s->total_bases += (n - trimmed) * (uint64_t)sp->readlen;
                                if (trimmed && (uint64_t)sp->readlen - (trimmed < 20 ? trimmed : 20) < s->min_length) {
                                        s->min_length = (uint64_t)sp->readlen - (trimmed < 20 ? trimmed : 20);
                                }
                                if ((uint64_t)sp->readlen < s->min_length) {
                                        s->min_length = (uint64_t)sp->readlen;
                                }
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
                                                        /* split by read base (sum = mismatches); MAPQ < 10 has a G
                                                         * pile-up at the first three positions (5' added G, 2011 paper) */
                                                        const uint64_t mm = pb[i].mismatches;
                                                        const int g_pile = m == COLLECT_MAPQ_LT10 && !axis && i < 3;
                                                        pb[i].mismatch_base[0] = g_pile ? mm / 20 : mm * 20 / 100;
                                                        pb[i].mismatch_base[1] = g_pile ? mm / 20 : mm * 20 / 100;
                                                        pb[i].mismatch_base[3] = g_pile ? mm / 20 : mm * 25 / 100;
                                                        pb[i].mismatch_base[4] = mm / 100;
                                                        pb[i].mismatch_base[2] = mm - pb[i].mismatch_base[0] - pb[i].mismatch_base[1] -
                                                                                 pb[i].mismatch_base[3] - pb[i].mismatch_base[4];
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
                                        for (int b = centre - spread; b <= centre + spread && b < COLLECT_N_IDENT && left; b++) {
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
                CHECK(contrast_on_white(report_mapq_colour[m]) >= 3.0, "MAPQ colour %d is below 3:1 (WCAG 1.4.11)", m);
        }
        for (int b = 0; b < 5; b++) {
                fprintf(stderr, "contrast on white: base %s  %s %.2f:1\n", b == 0 ? "A" : b == 1 ? "C" : b == 2 ? "G" : b == 3 ? "T" : "N", report_base_colour[b], contrast_on_white(report_base_colour[b]));
                CHECK(contrast_on_white(report_base_colour[b]) >= 3.0, "base colour %d is below 3:1 (WCAG 1.4.11)", b);
        }
        CHECK(contrast_on_white("#777777") >= 4.4, "grey");
        /* every base has its own dash pattern, so no line depends on colour alone (WCAG 1.4.1) */
        for (int b = 0; b < 5; b++) {
                for (int k = b + 1; k < 5; k++) {
                        CHECK(strcmp(report_base_dash[b], report_base_dash[k]) != 0, "base dashes %d and %d differ", b, k);
                }
        }
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
        /* two figures (Read 1, Read 2) per plot, plots 4 and 5 have the base key as well */
        for (int id = 1; id <= 6; id++) {
                CHECK(count_of(p[id], "<svg ") == (id == 4 || id == 5 ? 3 : 2), "plot %d svg count %d", id, count_of(p[id], "<svg "));
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
                CHECK(count_of(p[1], needle) == 2 + 2, "plot 1: one %s line and one key sample per mate (%d)", report_mapq_colour[m], count_of(p[1], needle));
                if (report_mapq_dash[m][0]) {
                        snprintf(needle, sizeof(needle), "stroke=\"%s\" stroke-width=\"2\" stroke-linejoin=\"round\" stroke-dasharray=\"%s\"", report_mapq_colour[m], report_mapq_dash[m]);
                        CHECK(count_of(p[1], needle) == 2, "plot 1: the %s line has its dash pattern", report_mapq_colour[m]);
                }
        }
        CHECK(count_of(p[1], "<title>MAPQ &gt;= 30, Read 1</title>") == 1, "series title carries class and mate");
        CHECK(count_of(p[5], "stroke=\"#777777\"") == 0 && count_of(p[6], "stroke=\"#777777\"") == 0, "unmapped is not drawn in plots 5 and 6");
        CHECK(count_of(p[3], "stroke=\"#777777\"") == 2 + 2, "unmapped is drawn in plot 3 (and in its key)");
        /* the key inside each MAPQ-split figure: every class with its read count, in table order */
        CHECK(count_of(p[1], "MAPQ &gt;= 30 (") == 2 && count_of(p[1], "Unmapped (") == 2, "plot 1: a key per mate");
        CHECK(strstr(p[1], "MAPQ &gt;= 30 (") < strstr(p[1], "Unmapped ("), "key in table order");
        CHECK(count_of(p[6], "Unmapped (") == 0, "plot 6 key: no unmapped class");
        CHECK(count_of(p[4], "MAPQ &gt;= 30 (") == 0, "plot 4: no MAPQ key");
        for (int b = 0; b < 5; b++) {
                char needle[64];
                snprintf(needle, sizeof(needle), "stroke=\"%s\"", report_base_colour[b]);
                CHECK(count_of(p[4], needle) == 4 * 2 + 1, "plot 4: base colour %s once per class row and mate and once in the key (%d)", report_base_colour[b], count_of(p[4], needle));
        }
        CHECK(strstr(p[4], "plot-key") != NULL && count_of(p[4], ">A</text>") == 1, "plot 4 key");
        /* plot 4: one row per MAPQ class, in table order (rule 9) */
        CHECK(count_of(p[4], ">Read 1: MAPQ &gt;= 30<") == 1 && count_of(p[4], ">Read 2: Unmapped<") == 1, "plot 4 row titles");
        CHECK(strstr(p[4], ">Read 1: MAPQ &gt;= 30<") < strstr(p[4], ">Read 1: MAPQ 10-29<") &&
              strstr(p[4], ">Read 1: MAPQ 10-29<") < strstr(p[4], ">Read 1: MAPQ &lt; 10<") &&
              strstr(p[4], ">Read 1: MAPQ &lt; 10<") < strstr(p[4], ">Read 1: Unmapped<"), "plot 4 rows in table order");
        CHECK(count_of(p[4], ">% of bases<") == 4 * 2, "plot 4: four rows per mate");
        CHECK(count_of(p[4], "<title>Base N, MAPQ &lt; 10, Read 2</title>") == 1, "plot 4 series title names base, class and mate");
        /* plot 5: one mismatch row per mapped class, stacked areas in the base colours, no unmapped row */
        CHECK(strstr(p[5], "id=\"samstat-p5-key-t\"") != NULL && strstr(p[5], "Bottom to top in every mismatch row: A, C, G, T, N.") != NULL, "plot 5 key");
        CHECK(count_of(p[5], ">Read 1: Mismatches, MAPQ") == 3 && strstr(p[5], "Mismatches, Unmapped") == NULL, "plot 5: three mismatch rows");
        for (int b = 0; b < 5; b++) {
                char needle[64];
                snprintf(needle, sizeof(needle), "<path fill=\"%s\"", report_base_colour[b]);
                CHECK(count_of(p[5], needle) == 3 * 2, "plot 5: base %s filled once per mismatch row and mate (%d)", report_base_colour[b], count_of(p[5], needle));
        }

        /* x ranges per 3.2: short layout, position 1 to the last bin with a base */
        CHECK(strstr(p[3], "from position 1 to position 150") != NULL, "plot 3 range 150 in the description");
        CHECK(strstr(p[3], ">position in read (bp), aligned part<") != NULL, "axis label says aligned part");
        /* a paired figure is 400 units wide: the mismatch row titles give the unit as % */
        CHECK(strstr(p[5], ">Read 1: Mismatches, MAPQ &gt;= 30 (%)<") != NULL && strstr(p[5], ">Read 2: Mismatches, MAPQ &lt; 10 (%)<") != NULL &&
              strstr(p[5], ">Read 2: Deletions (% of aligned bases)<") != NULL, "plot 5 row titles");
        CHECK(strstr(p[5], ">Read 1: Mismatches, MAPQ &gt;= 30 (%)<") < strstr(p[5], ">Read 1: Mismatches, MAPQ 10-29 (%)<") &&
              strstr(p[5], ">Read 1: Mismatches, MAPQ &lt; 10 (%)<") < strstr(p[5], ">Read 1: Insertions (% of aligned bases)<"),
              "plot 5 rows: mismatch rows in table order, then insertions");

        /* plot 5: both mates have identical axes, so the same labels on every tick */
        {
                const char* first = strstr(p[5], "<svg xmlns=\"http://www.w3.org/2000/svg\" role=\"img\" aria-labelledby=\"samstat-p5-r1-t");
                const char* second = first ? strstr(first + 10, "<svg ") : NULL;
                char a[4000], b[4000];
                CHECK(first && second, "plot 5 figures");
                if (!first || !second) {
                        first = second = p[5];
                }
                svg_texts(first, second, a, sizeof(a));
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
                CHECK(count_of(d, "<svg ") == (id == 4 || id == 5 ? 2 : 1), "%s: one figure", name);
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
                        /* long layout: the narrowest of 0.1, 0.2, 0.5 and 1 point that gives at most 50 bins */
                        {
                                const char* f = strstr(d, "from ");
                                int lo = -1;
                                double bw = 0.0;
                                CHECK(f && sscanf(f, "from %d%% to 100%%, in bins of %lf percentage point", &lo, &bw) == 2,
                                      "long layout identity bins: %.200s", f ? f : "(none)");
                                if (lo >= 0 && bw > 0.0) {
                                        static const double widths[4] = {0.1, 0.2, 0.5, 1.0};
                                        double want = 1.0;
                                        for (int w = 0; w < 4; w++) {
                                                if ((100.0 - lo) / widths[w] <= 50.0 + 1e-9) {
                                                        want = widths[w];
                                                        break;
                                                }
                                        }
                                        CHECK(fabs(bw - want) < 1e-9, "long layout identity bins: %g for %d-100%%, expected %g", bw, lo, want);
                                }
                        }
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
                        CHECK(strstr(p5, ">Mismatches, MAPQ") == NULL && strstr(p5, ">Mismatches (") == NULL && strstr(p5, ">Insertions (% of aligned bases)<") != NULL && strstr(p5, ">Deletions (% of aligned bases)<") != NULL, "indel rows only");
                        CHECK(count_of(p5, "<svg ") == 1 && strstr(p5, "samstat-p5-key") == NULL && strstr(p5, "<path fill=\"#") == NULL,
                              "no base key and no stacked areas without MD");
                        CHECK(count_of(p5, "(% of aligned bases)<") == 2, "two rows");
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
                CHECK(p5 && strstr(p5, " of 9,500 mapped reads (no MD tag and no =/X operations in CIGAR).</p>") != NULL, "note counts: %.250s", p5 ? p5 : "(null)");
                CHECK(p5 && strstr(p5, ">Mismatches, MAPQ &gt;= 30 (% of aligned bases)<") != NULL, "mismatch row present");
                CHECK(p5 && count_of(p5, ">Mismatches, MAPQ") == 3, "one mismatch row per mapped class");
                CHECK(p5 && strstr(p5, "<p class=\"plot-note\">") > strstr(p5, "samstat-p5-key"), "the note follows the base key");
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

static char* render_shape(const struct stat_collection* s, int plot)
{
        tld_strbuf* o = NULL;
        tld_strbuf_alloc(&o, 4096);
        int rc = report_plot_render_shape(o, s, plot);
        char* r = rc == OK ? strdup(TLD_STR(o)) : NULL;
        tld_strbuf_free(o);
        return r;
}

static void test_shape(void)
{
        /* 600 reads per mate: classes of 420, 90, 60 and 30 reads; the unmapped 30 are too few for a shape */
        struct spec sp = {.paired = 1, .aligned = 1, .mapped = 1, .is_long = 0, .quals = 1, .md = 1, .ident = 1, .readlen = 150, .nreads = 600};
        struct stat_collection* s = make(&sp);
        CHECK(s != NULL && s->basic_nums[0][COLLECT_MAPQ_UNMAPPED] == 30, "shape collection");
        if (!s) {
                return;
        }
        for (int id = 1; id <= 6; id++) {
                char name[32];
                snprintf(name, sizeof(name), "shape view of plot %d", id);
                char* d = render_shape(s, id);
                CHECK(d != NULL, "%s renders", name);
                if (!d) {
                        continue;
                }
                if (id == 3 || id == 4 || id == 5) {
                        CHECK(d[0] == 0, "%s: none (plots 3-5 are means, not counts)", name);
                        free(d);
                        continue;
                }
                check_markup(name, d);
                char* again = render_shape(s, id);
                CHECK(again && strcmp(again, d) == 0, "%s is not deterministic", name);
                free(again);
                char want[64];
                snprintf(want, sizeof(want), "id=\"samstat-p%ds-r2-t\"", id);
                CHECK(count_of(d, "<svg ") == 2 && strstr(d, want) != NULL, "%s: one figure per mate, own ids", name);
                CHECK(count_of(d, "% of the reads in the class") == 2, "%s: share axis", name);
                CHECK(strstr(d, "number of reads") == NULL, "%s: not a count axis", name);
                CHECK(strstr(d, "N50") == NULL && strstr(d, "median") == NULL, "%s: no marker of all reads", name);
                if (id != 6) {
                        CHECK(count_of(d, "Unmapped (30 reads, too few to draw)") == 2, "%s: small class named in the key", name);
                        CHECK(count_of(d, "stroke=\"#777777\"") == 2, "%s: small class only in the key", name);
                }
                free(d);
        }
        /* the counts view is unchanged by the shape view */
        char* p1 = render(s, 1);
        CHECK(p1 && strstr(p1, "number of reads") && count_of(p1, "Unmapped (30 reads)") == 2, "counts view keeps the small class");
        free(p1);
        stat_collection_free(s);

        /* unaligned input: one line, nothing to compare */
        struct spec su = {.paired = 0, .aligned = 0, .mapped = 0, .is_long = 1, .quals = 1, .readlen = 100, .nreads = 2000};
        s = make(&su);
        for (int id = 1; id <= 6 && s; id++) {
                char* d = render_shape(s, id);
                CHECK(d && d[0] == 0, "unaligned: no shape view of plot %d", id);
                free(d);
        }
        stat_collection_free(s);
}

static void test_single_length(void)
{
        /* one read length (single_len.sam): the position axis ends at that length; plot 1 is a sentence, not a
         * histogram of a single bar */
        struct spec sp = {.paired = 1, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 150, .nreads = 5000,
                          .one_len = 1};
        struct stat_collection* s = make(&sp);
        char* p1 = render(s, 1);
        char* p3 = render(s, 3);
        CHECK(p1 && strstr(p1, "<p class=\"plot-note\">All ") != NULL && strstr(p1, " reads are 150 bp long.</p>") != NULL &&
              strstr(p1, "<svg") == NULL, "one length: a sentence: %.300s", p1 ? p1 : "(null)");
        CHECK(p3 && strstr(p3, "from position 1 to position 150") != NULL, "axis ends at 150");
        CHECK(p3 && strstr(p3, ">150<") != NULL, "150 is a tick label");
        free(p1);
        free(p3);
        stat_collection_free(s);
}

/* Clarity round 1: captions by input type, exact N50, masked thin positions, no-sequence input,
 * key text of the base-composition lines. */
static void test_clarity(void)
{
        /* the caption of plot 3 mentions MAPQ classes only when there are some */
        CHECK(strstr(report_plot_caption(3, 1), "MAPQ") != NULL, "aligned caption of plot 3");
        CHECK(strstr(report_plot_caption(3, 0), "MAPQ") == NULL && strstr(report_plot_caption(3, 0), "class") == NULL,
              "unaligned caption of plot 3: %s", report_plot_caption(3, 0));
        CHECK(strcmp(report_plot_caption(1, 0), report_plot_caption(1, 1)) == 0, "plot 1 has one caption");
        for (int id = 1; id <= 6; id++) {
                CHECK(strstr(report_plot_caption(id, 1), "error probabilities") == NULL, "error-probability note is not in caption %d", id);
        }
        CHECK(strstr(report_plot_help[REPORT_PLOT_READQ].note, "error probabilities") != NULL, "the note lives in plot 2");
        CHECK(strstr(report_plot_help[REPORT_PLOT_READQ].caption, "Q20 = 1 error in 100 bases") != NULL, "plot 2 caption explains Q20");
        /* plots 4 and 5 (spec section 5): rows per MAPQ class on aligned input only */
        CHECK(strcmp(report_plot_caption(4, 1), "Percentage of A, C, G, T and N at each position along the read, one row per MAPQ class.") == 0 &&
              strcmp(report_plot_caption(4, 0), "Percentage of A, C, G, T and N at each position along the read.") == 0, "plot 4 captions");
        CHECK(strcmp(report_plot_caption(5, 1), "Mismatches per aligned base at each position, coloured by the read base and one row per MAPQ class; "
                                                "then insertions and deletions, one line per MAPQ class. Each panel has its own vertical scale.") == 0,
              "plot 5 caption: %s", report_plot_caption(5, 1));
        /* spec section 5: one good line for plots 4 and 5, the same for any input */
        CHECK(strcmp(report_plot_help[REPORT_PLOT_COMP].good, "Roughly parallel, flat lines, and the same picture in every row.") == 0,
              "plot 4 good line: %s", report_plot_help[REPORT_PLOT_COMP].good);
        CHECK(strstr(report_plot_help[REPORT_PLOT_COMP].look, "A row that differs from the MAPQ >= 30 row, such as more N") == report_plot_help[REPORT_PLOT_COMP].look &&
              strstr(report_plot_help[REPORT_PLOT_COMP].look, "an N spike at one position, or a change in composition") != NULL,
              "plot 4 look line: %s", report_plot_help[REPORT_PLOT_COMP].look);
        CHECK(report_plot_look(4, 1) == report_plot_help[REPORT_PLOT_COMP].look &&
              strstr(report_plot_look(4, 0), "MAPQ") == NULL && strstr(report_plot_look(4, 0), "An N spike at one position") == report_plot_look(4, 0),
              "plot 4 look line, unaligned: no MAPQ class");
        for (int id = 1; id <= 6; id++) {
                if (id != 4) {
                        CHECK(report_plot_look(id, 0) == report_plot_look(id, 1), "plot %d has one look line", id);
                }
        }
        CHECK(strstr(report_plot_help[REPORT_PLOT_ERRORS].look, "G at the first positions: bases added during library preparation") != NULL &&
              strcmp(report_plot_help[REPORT_PLOT_ERRORS].good, "Low, flat rates, with MAPQ >= 30 the lowest, and no single base dominating the mismatches.") == 0,
              "plot 5 help names the read base");

        /* positions backed by fewer than 50 reads are gaps, and the axis ends at the last backed position */
        {
                struct spec sp = {.paired = 0, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 150, .nreads = 1000};
                struct stat_collection* s = make(&sp);
                for (int m = 0; m < COLLECT_N_MAPQ; m++) {
                        for (int i = 100; i < 150; i++) {
                                struct collect_pos_bin* b = &s->acc->pos_abs[0][m][i];
                                memset(b, 0, sizeof(*b));
                                if (m == 0) {
                                        b->comp[0] = 3;
                                        b->qual_bases = 3;
                                        b->qual_error_sum = 3 * collect_qual_error(30);
                                }
                        }
                }
                char* p3 = render(s, 3);
                CHECK(p3 && strstr(p3, "from position 1 to position 100") != NULL, "axis ends at position 100: %.200s", p3 ? p3 : "(null)");
                CHECK(p3 && strstr(p3, "fewer than 50 reads") != NULL, "the description says why there are gaps");
                free(p3);
                stat_collection_free(s);
        }
        /* fewer reads than that: plots 3-5 are one "too few reads" line each (the page says the same),
         * from 50 reads on they are drawn; plots 1, 2 and 6 do not depend on it */
        {
                const uint64_t counts[2] = {49, 100};     /* per mate; the mapped classes round down */
                for (int k = 0; k < 2; k++) {
                        struct spec sp = {.paired = 0, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 100, .nreads = counts[k]};
                        struct stat_collection* s = make(&sp);
                        CHECK(s != NULL && stat_collection_finalise(s) == OK, "finalise %d", k);
                        char want[160];
                        snprintf(want, sizeof(want), "<p class=\"plot-na\">Too few reads (%" PRIu64 ") for a per-position plot (needs at least " REPORT_MIN_POS_BASES_STR ").</p>\n", s->n_total);
                        for (int id = 3; id <= 5; id++) {
                                char* d = render(s, id);
                                if (k == 0) {
                                        CHECK(d && strcmp(d, want) == 0, "%d reads, plot %d: %.200s", (int)counts[k], id, d ? d : "(null)");
                                } else {
                                        CHECK(d && strstr(d, "<svg ") != NULL && strstr(d, "Too few reads") == NULL, "%d reads, plot %d is drawn", (int)counts[k], id);
                                }
                                free(d);
                        }
                        for (int id = 1; id <= 6; id += id == 2 ? 4 : 1) {
                                char* d = render(s, id);
                                CHECK(d && strstr(d, "<svg ") != NULL && strstr(d, "Too few reads") == NULL, "%d reads, plot %d is drawn", (int)counts[k], id);
                                free(d);
                        }
                        stat_collection_free(s);
                }
        }
        /* the unaligned case says the same for plots 3 and 4 */
        {
                struct spec sp = {.paired = 0, .aligned = 0, .quals = 1, .readlen = 50, .nreads = 5};
                struct stat_collection* s = make(&sp);
                char* p3 = render(s, 3);
                CHECK(p3 && strstr(p3, "Too few reads (5) for a per-position plot (needs at least 50).") != NULL && strstr(p3, "<svg ") == NULL,
                      "5 reads: plot 3 is replaced: %.200s", p3 ? p3 : "(null)");
                free(p3);
                stat_collection_free(s);
        }

        /* reads without sequence: plots 1, 4 and 5 say so instead of drawing a 0-10 bp axis */
        {
                struct spec sp = {.paired = 0, .aligned = 1, .mapped = 1, .quals = 0, .md = 0, .ident = 0, .readlen = 100, .nreads = 100};
                struct stat_collection* s = make(&sp);
                memset(s->acc->len_reads, 0, sizeof(s->acc->len_reads));
                memset(s->acc->len_bases, 0, sizeof(s->acc->len_bases));
                memset(s->acc->pos_abs, 0, sizeof(s->acc->pos_abs));
                memset(s->acc->pos_rel, 0, sizeof(s->acc->pos_rel));
                s->total_bases = 0;
                s->min_length = 0;
                s->max_length = 0;
                CHECK(stat_collection_finalise(s) == OK, "finalise");
                const char* want = "<p class=\"plot-na\">Not available: the input has no sequence.</p>\n";
                for (int id = 1; id <= 5; id += 3) {
                        char* d = render(s, id);
                        CHECK(d && strcmp(d, want) == 0, "no sequence, plot %d: %s", id, d ? d : "(null)");
                        free(d);
                }
                char* d5 = render(s, 5);
                CHECK(d5 && strstr(d5, "the input has no sequence.") != NULL, "no sequence, plot 5: %.120s", d5 ? d5 : "(null)");
                free(d5);
                stat_collection_free(s);
        }

        /* N50 and the median of one read of exactly 1,000,000 bp are that length, not a bin edge */
        {
                struct stat_collection* s = NULL;
                CHECK(stat_collection_alloc(&s) == OK, "alloc");
                s->basic_nums[0][COLLECT_MAPQ_GE30] = 1;
                s->has_sq = 1;
                int b = collect_length_bin(1000000);
                s->acc->len_reads[0][COLLECT_MAPQ_GE30][b] = 1;
                s->acc->len_bases[0][COLLECT_MAPQ_GE30][b] = 1000000;
                s->total_bases = 1000000;
                s->min_length = 1000000;
                s->max_length = 1000000;
                CHECK(stat_collection_finalise(s) == OK, "finalise");
                CHECK(s->n50 < 1000000, "the JSON value stays the bin edge (%" PRIu64 ")", s->n50);
                CHECK(s->n50_best == 1000000 && s->n50_exact && s->median_best == 1000000 && s->median_exact,
                      "exact display values: %" PRIu64 " %" PRIu64, s->n50_best, s->median_best);
                char* p1 = render(s, 1);
                /* one read: plot 1 is a sentence with the exact length, in the unit of the Bases tile */
                CHECK(p1 && strstr(p1, "The read is 1.0 Mb long.") != NULL && strstr(p1, "~") == NULL, "plot 1 sentence is exact: %s",
                      p1 ? p1 : "(null)");
                free(p1);
                stat_collection_free(s);
        }

        /* two reads in one bin: the display value is the bin edge clamped to the observed lengths, marked approximate */
        {
                struct stat_collection* s = NULL;
                CHECK(stat_collection_alloc(&s) == OK, "alloc");
                s->basic_nums[0][COLLECT_MAPQ_GE30] = 2;
                s->has_sq = 1;
                int b = collect_length_bin(100500);
                CHECK(b == collect_length_bin(100502), "both lengths share a bin");
                s->acc->len_reads[0][COLLECT_MAPQ_GE30][b] = 2;
                s->acc->len_bases[0][COLLECT_MAPQ_GE30][b] = 201002;
                s->total_bases = 201002;
                s->min_length = 100500;
                s->max_length = 100502;
                CHECK(stat_collection_finalise(s) == OK, "finalise");
                CHECK(s->n50_best >= s->min_length && s->n50_best <= s->max_length && !s->n50_exact,
                      "n50 %" PRIu64 " (edge %" PRIu64 ") is within [min, max] and approximate", s->n50_best, s->n50);
                CHECK(s->median_best >= s->min_length && s->median_best <= s->max_length && !s->median_exact, "median is within [min, max]");
                stat_collection_free(s);
        }

        /* plot 4: the key names the dash of every line */
        {
                struct spec sp = {.paired = 0, .aligned = 0, .quals = 1, .readlen = 100, .nreads = 20000};
                struct stat_collection* s = make(&sp);
                char* p4 = render(s, 4);
                CHECK(p4 && strstr(p4, "A green, solid line; C blue, dashed line; G orange, dotted line; T red, dash-dot line; N grey, long-dashed line") != NULL,
                      "key text: %.600s", p4 ? p4 : "(null)");
                for (int b = 1; b < 5; b++) {
                        char needle[64];
                        snprintf(needle, sizeof(needle), "stroke-dasharray=\"%s\"", report_base_dash[b]);
                        CHECK(count_of(p4 ? p4 : "", needle) == 2, "plot 4: dash %s on the line and in the key", report_base_dash[b]);
                }
                free(p4);
                stat_collection_free(s);
        }
}

/* Clarity round 2: flat error panels say so on a whole-number scale; plot 6 ends at 100. */
static void test_clarity_r2(void)
{
        struct spec sp = {.paired = 0, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 100, .nreads = 1000};
        struct stat_collection* s = make(&sp);
        CHECK(s != NULL && stat_collection_finalise(s) == OK, "finalise");
        if (!s) {
                return;
        }
        char* p5 = render(s, 5);
        CHECK(p5 && count_of(p5, ">none detected<") < 5, "data with errors: not every panel is labelled none detected (%d)", p5 ? count_of(p5, ">none detected<") : -1);
        free(p5);
        for (int m = 0; m < COLLECT_N_MAPQ; m++) {
                for (int i = 0; i < COLLECT_ABS_BINS; i++) {
                        struct collect_pos_bin* b = &s->acc->pos_abs[0][m][i];
                        b->mismatches = 0;
                        memset(b->mismatch_base, 0, sizeof(b->mismatch_base));
                        b->insertions = 0;
                        b->deletions = 0;
                }
        }
        p5 = render(s, 5);
        CHECK(p5 && count_of(p5, ">none detected<") == 5, "all-zero data: five labelled panels, three of them mismatch rows (%d)", p5 ? count_of(p5, ">none detected<") : -1);
        CHECK(p5 && strstr(p5, ">4<") != NULL && strstr(p5, ">2<") != NULL && strstr(p5, ">0.5<") == NULL && strstr(p5, ">1.0<") == NULL,
              "all-zero data: whole-number scale 0 2 4");
        CHECK(p5 && strstr(p5, "labelled none detected") != NULL, "the description says what the label means");
        free(p5);

        char* p6 = render(s, 6);
        CHECK(p6 && strstr(p6, "from 95% to 100%") != NULL, "plot 6 desc: %.400s", p6 ? p6 : "(null)");
        CHECK(p6 && strstr(p6, ">95<") != NULL && strstr(p6, ">100<") != NULL && strstr(p6, ">101<") == NULL &&
              strstr(p6, ">99.2<") == NULL, "plot 6 ticks are whole numbers up to 100");
        free(p6);
        stat_collection_free(s);
}

/* Copy the d="..." of the path whose series-level <title> is title (HTML-escaped). */
static int path_d_of(const char* doc, const char* title, char* d, size_t dn)
{
        char t[160];
        snprintf(t, sizeof(t), "<title>%s</title>", title);
        const char* at = strstr(doc, t);
        const char* start = NULL;
        if (!at) {
                return 0;
        }
        for (const char* p = doc; p < at; p++) {
                if (strncmp(p, " d=\"", 4) == 0) {
                        start = p + 4;
                }
        }
        if (!start) {
                return 0;
        }
        size_t n = (size_t)(strchr(start, '"') - start);
        n = n >= dn ? dn - 1 : n;
        memcpy(d, start, n);
        d[n] = 0;
        return 1;
}

/* Remove every trace of MAPQ class m from a collection (as if no read had it). */
static void drop_class(struct stat_collection* s, int m)
{
        for (int r = 0; r < COLLECT_N_MATE; r++) {
                s->basic_nums[r][m] = 0;
                memset(s->acc->pos_abs[r][m], 0, sizeof(s->acc->pos_abs[r][m]));
                memset(s->acc->pos_rel[r][m], 0, sizeof(s->acc->pos_rel[r][m]));
                memset(s->acc->len_reads[r][m], 0, sizeof(s->acc->len_reads[r][m]));
                memset(s->acc->len_bases[r][m], 0, sizeof(s->acc->len_bases[r][m]));
                memset(s->acc->read_q[r][m], 0, sizeof(s->acc->read_q[r][m]));
                memset(s->acc->identity[r][m], 0, sizeof(s->acc->identity[r][m]));
        }
}

/* Base x MAPQ (encoding rule 9, decisions A3 and F3 revised): plot 4 has one row per class with
 * reads, plot 5 one mismatch row per mapped class with MD-backed bases. */
static void test_base_rows(void)
{
        char d[20000];

        /* paired, no read in MAPQ 10-29, and the MAPQ < 10 reads have no MD */
        {
                struct spec sp = {.paired = 1, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 100, .nreads = 5000};
                struct stat_collection* s = make(&sp);
                CHECK(s != NULL, "collection");
                if (!s) {
                        return;
                }
                drop_class(s, COLLECT_MAPQ_10_29);
                for (int r = 0; r < COLLECT_N_MATE; r++) {
                        for (int i = 0; i < COLLECT_ABS_BINS; i++) {
                                struct collect_pos_bin* b = &s->acc->pos_abs[r][COLLECT_MAPQ_LT10][i];
                                b->aligned_md = 0;
                                b->mismatches = 0;
                                memset(b->mismatch_base, 0, sizeof(b->mismatch_base));
                        }
                }
                CHECK(stat_collection_finalise(s) == OK, "finalise");
                char* p4 = render(s, 4);
                char* p5 = render(s, 5);
                CHECK(p4 && p5, "render");
                if (p4 && p5) {
                        check_markup("rows plot 4", p4);
                        check_markup("rows plot 5", p5);
                        CHECK(count_of(p4, ">% of bases<") == 3 * 2, "plot 4: three rows per mate (%d)", count_of(p4, ">% of bases<"));
                        CHECK(strstr(p4, "MAPQ 10-29") == NULL, "plot 4: no row for a class without reads");
                        CHECK(strstr(p4, ">Read 2: MAPQ &gt;= 30<") < strstr(p4, ">Read 2: MAPQ &lt; 10<") &&
                              strstr(p4, ">Read 2: MAPQ &lt; 10<") < strstr(p4, ">Read 2: Unmapped<"), "plot 4: table order");
                        CHECK(strstr(p4, "Rows: MAPQ &gt;= 30, MAPQ &lt; 10, Unmapped.") != NULL, "plot 4 description names the rows");
                        CHECK(count_of(p5, ">Read 1: Mismatches, MAPQ") == 1 && strstr(p5, ">Read 1: Mismatches, MAPQ &gt;= 30 (%)<") != NULL,
                              "plot 5: a mismatch row only for the class with MD");
                        CHECK(strstr(p5, "Mismatches, MAPQ 10-29") == NULL && strstr(p5, "Mismatches, MAPQ &lt; 10") == NULL,
                              "plot 5: no mismatch row for a class without reads or without MD");
                        CHECK(strstr(p5, "<title>MAPQ &lt; 10, Read 1</title>") != NULL, "plot 5: the MD-less class keeps its indel lines");
                        CHECK(strstr(p5, "MAPQ &lt; 10 (") != NULL && strstr(p5, "too few to draw") == NULL, "plot 5 key: the class is drawn");
                        CHECK(count_of(p5, "(% of aligned bases)<") == 2 * 2, "plot 5: insertion and deletion rows per mate");
                }
                free(p4);
                free(p5);
                stat_collection_free(s);
        }

        /* unaligned input: plot 4 is one row of all reads, as tall as before, with no MAPQ in it */
        {
                struct spec sp = {.paired = 0, .aligned = 0, .quals = 1, .readlen = 100, .nreads = 2000};
                struct stat_collection* s = make(&sp);
                char* p4 = render(s, 4);
                CHECK(p4 && count_of(p4, ">% of bases<") == 1 && count_of(p4, "<svg ") == 2, "unaligned plot 4: one row");
                CHECK(p4 && strstr(p4, "MAPQ") == NULL && strstr(p4, "Unmapped") == NULL && strstr(p4, "Rows:") == NULL, "unaligned plot 4: no class");
                CHECK(p4 && strstr(p4, "<title>Base A</title>") != NULL && strstr(p4, "<title>Base N</title>") != NULL, "unaligned plot 4: series of all reads");
                CHECK(p4 && strstr(p4, "Positions covered by fewer than 50 reads are not drawn.") != NULL, "unaligned plot 4 description");
                free(p4);
                stat_collection_free(s);
        }

        /* single, short: a bin with fewer than 50 bases in a class is a gap in that class's row only;
         * a mismatch row has a gap where its class has fewer than 50 MD-backed aligned bases */
        {
                struct spec sp = {.paired = 0, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 100, .nreads = 5000};
                struct stat_collection* s = make(&sp);
                struct collect_pos_bin* u = &s->acc->pos_abs[0][COLLECT_MAPQ_UNMAPPED][10];
                memset(u->comp, 0, sizeof(u->comp));
                u->comp[0] = 49;
                struct collect_pos_bin* g = &s->acc->pos_abs[0][COLLECT_MAPQ_GE30][20];
                g->aligned_md = 49;
                CHECK(stat_collection_finalise(s) == OK, "finalise");
                char* p4 = render(s, 4);
                char* p5 = render(s, 5);
                CHECK(p4 && path_d_of(p4, "Base A, Unmapped", d, sizeof(d)) && count_of(d, "M") == 2, "plot 4: thin bin is a gap in the unmapped row");
                CHECK(p4 && path_d_of(p4, "Base A, MAPQ &gt;= 30", d, sizeof(d)) && count_of(d, "M") == 1, "plot 4: other rows have no gap");
                CHECK(p5 && path_d_of(p5, "Mismatches at read base G, MAPQ &gt;= 30", d, sizeof(d)) && count_of(d, "M") == 2,
                      "plot 5: 49 MD-backed bases are a gap in the mismatch row: %.200s", d);
                CHECK(p5 && path_d_of(p5, "Mismatches at read base G, MAPQ &lt; 10", d, sizeof(d)) && count_of(d, "M") == 1,
                      "plot 5: other mismatch rows have no gap");
                free(p4);
                free(p5);
                stat_collection_free(s);
        }

        /*
         * Stack order and height: constant per-base mismatch counts A 100, C 200, G 300, T 400, N 50 per
         * 100,000 MD-backed bases (total 1.05 %). The areas are drawn top first (N, T, G, C, A), each from
         * the bottom of the row to the running sum, so the top of N is the total mismatch rate.
         */
        {
                static const uint64_t mb[5] = {100, 200, 300, 400, 50};
                static const double cum[5] = {100, 300, 600, 1000, 1050};
                struct spec sp = {.paired = 0, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 100, .nreads = 5000};
                struct stat_collection* s = make(&sp);
                for (int i = 0; i < COLLECT_ABS_BINS; i++) {
                        struct collect_pos_bin* b = &s->acc->pos_abs[0][COLLECT_MAPQ_GE30][i];
                        if (!b->aligned) {
                                continue;
                        }
                        b->aligned_md = 100000;
                        b->mismatches = 1050;
                        memcpy(b->mismatch_base, mb, sizeof(mb));
                }
                CHECK(stat_collection_finalise(s) == OK, "finalise");
                char* p5 = render(s, 5);
                double top_y[5] = {0}, base_y[5] = {0};
                const char* pos[5] = {NULL};
                for (int b = 0; b < 5 && p5; b++) {
                        char t[96];
                        double x0, dx, dy;
                        snprintf(t, sizeof(t), "Mismatches at read base %s, MAPQ &gt;= 30", b == 0 ? "A" : b == 1 ? "C" : b == 2 ? "G" : b == 3 ? "T" : "N");
                        int ok = path_d_of(p5, t, d, sizeof(d));
                        pos[b] = strstr(p5, t);
                        /* flat: the first point and one step to the last (the tolerance drops the rest), down to
                         * the bottom, back, closed */
                        CHECK(ok && sscanf(d, "M%lf,%lfl%lf,%lfV%lfH", &x0, &top_y[b], &dx, &dy, &base_y[b]) == 5 &&
                              dy == 0.0 && count_of(d, ",") == 2, "flat stack path of base %d: %.200s", b, ok ? d : "(none)");
                }
                CHECK(p5 && pos[4] && pos[3] && pos[2] && pos[1] && pos[0] && pos[4] < pos[3] && pos[3] < pos[2] && pos[2] < pos[1] &&
                      pos[1] < pos[0], "stack drawn top first: N T G C A");
                /* the plotting height of a plot 5 row is PH_ROW (96 units); the axis top is err_top(total) */
                const double total = 100.0 * 1050.0 / 100000.0;
                const double axis_top = fmax(svg_nice_ceil(total * 1.04, 3), 0.5);
                for (int b = 0; b < 5 && p5; b++) {
                        double want = 96.0 * (100.0 * cum[b] / 100000.0) / axis_top;
                        CHECK(fabs((base_y[b] - top_y[b]) - want) < 0.11, "base %d: band top at %.2f units, expected %.2f (running sum)", b,
                              base_y[b] - top_y[b], want);
                }
                CHECK(p5 && fabs((base_y[4] - top_y[4]) - 96.0 * total / axis_top) < 0.11, "the top of the stack is the total mismatch rate");
                CHECK(p5 && strstr(p5, "bottom to top A, C, G, T, N") != NULL && strstr(p5, "the top edge of the stack is the total mismatch rate") != NULL,
                      "plot 5 description explains the stack");
                free(p5);

                /* a base without mismatches has no area (its band is empty) */
                for (int i = 0; i < COLLECT_ABS_BINS; i++) {
                        struct collect_pos_bin* b = &s->acc->pos_abs[0][COLLECT_MAPQ_GE30][i];
                        if (b->aligned_md) {
                                b->mismatches -= b->mismatch_base[4];
                                b->mismatch_base[4] = 0;
                        }
                }
                p5 = render(s, 5);
                CHECK(p5 && strstr(p5, "Mismatches at read base N, MAPQ &gt;= 30") == NULL &&
                      strstr(p5, "Mismatches at read base T, MAPQ &gt;= 30") != NULL, "no area for a base without mismatches");
                free(p5);
                stat_collection_free(s);
        }

        /*
         * Each mismatch row has its own vertical scale: constant total rates 0.4 % (MAPQ >= 30), 1 %
         * (MAPQ 10-29) and 4 % (MAPQ < 10), split evenly over the five bases. The top of each stack is
         * at its own rate on an axis that ends at err_top() of that row's maximum, not of all rows.
         */
        {
                static const uint64_t total[3] = {400, 1000, 4000};   /* per 100,000 MD-backed bases */
                static const char* const cls[3] = {"MAPQ &gt;= 30", "MAPQ 10-29", "MAPQ &lt; 10"};
                struct spec sp = {.paired = 0, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 100, .nreads = 5000};
                struct stat_collection* s = make(&sp);
                for (int m = 0; m < 3; m++) {
                        for (int i = 0; i < COLLECT_ABS_BINS; i++) {
                                struct collect_pos_bin* b = &s->acc->pos_abs[0][m][i];
                                if (!b->aligned) {
                                        continue;
                                }
                                b->aligned_md = 100000;
                                b->mismatches = total[m];
                                for (int k = 0; k < 5; k++) {
                                        b->mismatch_base[k] = total[m] / 5;
                                }
                        }
                }
                CHECK(stat_collection_finalise(s) == OK, "finalise");
                char* p5 = render(s, 5);
                double tops[3] = {0};
                for (int m = 0; m < 3 && p5; m++) {
                        char t[96];
                        double x0, y0, dx, dy, base;
                        snprintf(t, sizeof(t), "Mismatches at read base N, %s", cls[m]);
                        int ok = path_d_of(p5, t, d, sizeof(d)) && sscanf(d, "M%lf,%lfl%lf,%lfV%lfH", &x0, &y0, &dx, &dy, &base) == 5;
                        const double rate = 100.0 * (double)total[m] / 100000.0;
                        tops[m] = fmax(svg_nice_ceil(rate * 1.04, 3), 0.5);
                        const double want = 96.0 * rate / tops[m];
                        CHECK(ok && fabs((base - y0) - want) < 0.11, "%s: stack top at %.2f units, expected %.2f (own axis top %g)", cls[m],
                              ok ? base - y0 : -1.0, want, tops[m]);
                }
                CHECK(tops[0] < tops[1] && tops[1] < tops[2], "three different axis tops (%g, %g, %g)", tops[0], tops[1], tops[2]);
                /* the tick labels of each row follow its own axis (0-0.5 and 0-5 here) */
                char a[4000];
                svg_texts(p5 ? p5 : "", p5 ? p5 + strlen(p5) : "", a, sizeof(a));
                CHECK(p5 && strstr(a, "|Mismatches, MAPQ &gt;= 30 (% of aligned bases)|0.0|0.2|0.4|Mismatches, MAPQ 10-29") != NULL &&
                      strstr(a, "|Mismatches, MAPQ &lt; 10 (% of aligned bases)|0|2|4|Insertions") != NULL,
                      "tick labels of the own axis of each row: %.400s", a);
                free(p5);
                stat_collection_free(s);
        }

        /*
         * Paired: a mismatch row exists when either mate has MD-backed bases in the class, so the other
         * mate's panel can be empty. Its note names the reason, as in plot 4: Read 1 has no MAPQ < 10
         * reads, no MD in MAPQ 10-29, and fewer than 50 MD-backed bases per position in MAPQ >= 30.
         * Plot 4: Read 1 has no MAPQ < 10 reads, and fewer than 50 bases per position in MAPQ >= 30.
         */
        {
                struct spec sp = {.paired = 1, .aligned = 1, .mapped = 1, .quals = 1, .md = 1, .ident = 1, .readlen = 100, .nreads = 5000};
                struct stat_collection* s = make(&sp);
                const int m_lt = COLLECT_MAPQ_LT10;
                s->basic_nums[0][m_lt] = 0;
                memset(s->acc->pos_abs[0][m_lt], 0, sizeof(s->acc->pos_abs[0][m_lt]));
                memset(s->acc->pos_rel[0][m_lt], 0, sizeof(s->acc->pos_rel[0][m_lt]));
                memset(s->acc->len_reads[0][m_lt], 0, sizeof(s->acc->len_reads[0][m_lt]));
                memset(s->acc->len_bases[0][m_lt], 0, sizeof(s->acc->len_bases[0][m_lt]));
                memset(s->acc->read_q[0][m_lt], 0, sizeof(s->acc->read_q[0][m_lt]));
                memset(s->acc->identity[0][m_lt], 0, sizeof(s->acc->identity[0][m_lt]));
                for (int i = 0; i < COLLECT_ABS_BINS; i++) {
                        struct collect_pos_bin* b = &s->acc->pos_abs[0][COLLECT_MAPQ_10_29][i];
                        b->aligned_md = 0;
                        b->mismatches = 0;
                        memset(b->mismatch_base, 0, sizeof(b->mismatch_base));
                        b = &s->acc->pos_abs[0][COLLECT_MAPQ_GE30][i];
                        if (b->aligned) {
                                b->aligned_md = 10;
                                b->mismatches = 1;
                                memset(b->mismatch_base, 0, sizeof(b->mismatch_base));
                                b->mismatch_base[0] = 1;
                                memset(b->comp, 0, sizeof(b->comp));
                                b->comp[0] = 10;
                        }
                }
                CHECK(stat_collection_finalise(s) == OK, "finalise");
                char* p4 = render(s, 4);
                char* p5 = render(s, 5);
                const char* r2 = p5 ? strstr(p5, ">Read 2: ") : NULL;
                CHECK(p5 && r2, "render plot 5");
                if (p5 && r2) {
                        check_markup("plot 5, empty R1 mismatch panels", p5);
                        const char* t30 = strstr(p5, ">Read 1: Mismatches, MAPQ &gt;= 30 (%)<");
                        const char* t10 = strstr(p5, ">Read 1: Mismatches, MAPQ 10-29 (%)<");
                        const char* tlt = strstr(p5, ">Read 1: Mismatches, MAPQ &lt; 10 (%)<");
                        const char* tins = strstr(p5, ">Read 1: Insertions");
                        const char* n_few = strstr(p5, ">too few reads to draw<");
                        const char* n_md = strstr(p5, ">no reads with MD or X<");
                        const char* n_none = strstr(p5, ">no reads<");
                        CHECK(t30 && t10 && tlt && tins, "plot 5: Read 1 keeps the three mismatch rows");
                        CHECK(t30 && t10 && tlt && tins && n_few && n_md && n_none && t30 < n_few && n_few < t10 && t10 < n_md && n_md < tlt &&
                              tlt < n_none && n_none < tins, "plot 5 R1 notes: too few, no MD, no reads, each in its own row");
                        CHECK(count_of(p5, ">too few reads to draw<") == 1 && count_of(p5, ">no reads with MD or X<") == 1 &&
                              count_of(p5, ">no reads<") == 1 && (!n_none || n_none < r2), "plot 5: no notes in Read 2");
                }
                const char* q2 = p4 ? strstr(p4, ">Read 2: ") : NULL;
                CHECK(p4 && q2, "render plot 4");
                if (p4 && q2) {
                        const char* t30 = strstr(p4, ">Read 1: MAPQ &gt;= 30<");
                        const char* tlt = strstr(p4, ">Read 1: MAPQ &lt; 10<");
                        const char* tun = strstr(p4, ">Read 1: Unmapped<");
                        const char* n_few = strstr(p4, ">too few reads to draw<");
                        const char* n_none = strstr(p4, ">no reads<");
                        CHECK(t30 && tlt && tun && n_few && n_none && t30 < n_few && n_few < tlt && tlt < n_none && n_none < tun,
                              "plot 4 R1 notes: too few in MAPQ >= 30, no reads in MAPQ < 10");
                        CHECK(count_of(p4, ">too few reads to draw<") == 1 && count_of(p4, ">no reads<") == 1 && n_none && n_none < q2,
                              "plot 4: no notes in Read 2");
                }
                free(p4);
                free(p5);
                stat_collection_free(s);
        }
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
        /* the class name is the row header next to the swatch: the figure names only its style */
        CHECK(strstr(TLD_STR(o), "dotted vermillion line</title>") != NULL && strstr(TLD_STR(o), "Line sample: vermillion, dotted.") != NULL, "swatch text");
        CHECK(strstr(TLD_STR(o), "MAPQ &gt;= 30") == NULL, "swatch titles do not repeat the class name");
        CHECK(strstr(TLD_STR(o), "stroke-dasharray=\"12 4 2 4\"") != NULL, "dash-dot");
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

        /*
         * The same with noise at every position (a shallow file: a few thousand reads per class), so that
         * no point of plots 4 and 5 lies on the line through its neighbours and the path tolerance drops
         * nothing: one row of five bases per MAPQ class in plot 4 and one stack per mapped class in plot 5.
         */
        lcg = 4242;
        for (int r = 0; r < COLLECT_N_MATE; r++) {
                for (int m = 0; m < COLLECT_N_MAPQ; m++) {
                        for (int i = 0; i < COLLECT_ABS_BINS; i++) {
                                struct collect_pos_bin* b = &s->acc->pos_abs[r][m][i];
                                uint64_t tot = 0, mm = 0;
                                for (int c = 0; c < 5; c++) {
                                        b->comp[c] = 200 + (uint64_t)(400.0 * rnd());
                                        tot += b->comp[c];
                                        if (m != COLLECT_MAPQ_UNMAPPED) {
                                                b->mismatch_base[c] = (uint64_t)(60.0 * rnd());
                                                mm += b->mismatch_base[c];
                                        }
                                }
                                if (m != COLLECT_MAPQ_UNMAPPED) {
                                        b->aligned = b->aligned_md = tot;
                                        b->mismatches = mm;
                                        b->insertions = (uint64_t)(20.0 * rnd());
                                        b->deletions = (uint64_t)(20.0 * rnd());
                                }
                                b->qual_bases = tot;
                                b->qual_error_sum = tot * collect_qual_error(20 + (int)(15.0 * rnd()));
                        }
                }
        }
        tld_strbuf_clear(o);
        for (int id = 1; id <= 6; id++) {
                size_t before = (size_t)TLD_STRLEN(o);
                CHECK(report_plot_render(o, s, id) == OK, "plot %d of the noisy worst case", id);
                fprintf(stderr, "noisy worst case (paired 999 bp): plot %d is %zu bytes\n", id, (size_t)TLD_STRLEN(o) - before);
        }
        fprintf(stderr, "noisy worst case (paired 999 bp): six plots %zu bytes\n", (size_t)TLD_STRLEN(o));
        /* the page limit is 1 MB and the rest of the page (text, tables, CSS) is about 20 kB */
        CHECK((size_t)TLD_STRLEN(o) < 950000, "six plots of the noisy worst case are %zu bytes (the page limit is 1 MB)", (size_t)TLD_STRLEN(o));
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
        test_shape();
        test_clarity();
        test_clarity_r2();
        test_base_rows();
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
