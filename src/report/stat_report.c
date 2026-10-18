#include "tldlite.h"

#include <inttypes.h>
#include <math.h>
#include <string.h>
#include "hts.h"
#include "../plot/report_plots.h"
#include "../plot/report_fmt.h"
#include "../param/param.h"
#include "../collect/collect.h"
#include "config.h"
#include "report_help.h"
#include "outfile.h"
#define  STAT_REPORT_IMPORT
#include "stat_report.h"
#include "json_report.h"
#include "mqc_report.h"

/*
 * The standard report page (docs/plan-v2.4/standard-plots.md section 1).
 *
 * One self-contained HTML file, no script, no external resource:
 *   header, summary strip, notes, MAPQ table (aligned input), plots 1-6, footer.
 * All numbers come from the finalised collection (stat_collection_finalise); nothing
 * is recomputed here except presentation (rounding, units). All wording that is not
 * layout lives in report_help.c, the figures in report_plots.c.
 */

#define MAX_NUM 64          /* longest formatted number */

static int page_head(tld_strbuf *o, const char *name);
static int page_header(tld_strbuf *o, const char *name, const char *format, struct stat_collection *s);
static int page_strip(tld_strbuf *o, struct stat_collection *s);
static int page_notes(tld_strbuf *o, struct stat_collection *s);
static int page_mapq_table(tld_strbuf *o, struct stat_collection *s);
static int page_plot(tld_strbuf *o, struct stat_collection *s, int plot_id);
static int page_footer(tld_strbuf *o, const char *name, struct stat_collection *s);

static int write_buffer(tld_strbuf *out_buffer, struct samstat_param *p, int id);

/* ------------------------------------------------------------------ */
/* small helpers                                                        */
/* ------------------------------------------------------------------ */

static int app(tld_strbuf *o, const char *text)
{
        return tld_append(o, (char *) text);
}

/* text with &, <, > and quotes escaped */
static int esc(tld_strbuf *o, const char *text)
{
        for(const char *c = text; *c; c++){
                switch(*c){
                case '&':
                        RUN(app(o, "&amp;"));
                        break;
                case '<':
                        RUN(app(o, "&lt;"));
                        break;
                case '>':
                        RUN(app(o, "&gt;"));
                        break;
                case '"':
                        RUN(app(o, "&quot;"));
                        break;
                case '\'':
                        RUN(app(o, "&#39;"));
                        break;
                default:
                        RUN(tld_append_char(o, *c));
                        break;
                }
        }
        return OK;
ERROR:
        return FAIL;
}

/* A template from report_help.c: every "%s" is replaced by the next argument (numbers
 * formatted by the caller, so they need no escaping), "%%" is a percent sign, the rest is
 * escaped text. Not printf: the templates are data. */
static int app_template(tld_strbuf *o, const char *t, const char *a1, const char *a2)
{
        const char *arg[2] = {a1, a2};
        int n_used = 0;
        for(const char *c = t; *c; c++){
                if(c[0] == '%' && c[1] == 's'){
                        if(n_used < 2 && arg[n_used]){
                                RUN(app(o, arg[n_used]));
                        }
                        n_used++;
                        c++;
                }else if(c[0] == '%' && c[1] == '%'){
                        RUN(app(o, "%"));
                        c++;
                }else{
                        char one[2] = {*c, 0};
                        RUN(esc(o, one));
                }
        }
        return OK;
ERROR:
        return FAIL;
}

/* Counts and bases: exact with thousands separators below 10,000, then one decimal and an SI
 * suffix (12.4 M reads, 12.4 Mb). unit is "" for reads and "b" for bases; bases below 10,000 get
 * " bp". */
static void fmt_si_unit(char *buf, size_t n, uint64_t v, const char *unit)
{
        static const char *const prefix[] = {"k", "M", "G", "T", "P", "E"};
        double x;
        int u = -1;
        char d[32];

        if(v < 10000){
                fmt_count(buf, n, v);
                if(unit[0]){
                        size_t l = strlen(buf);
                        snprintf(buf + l, n - l, " bp");
                }
                return;
        }
        x = (double) v / 1000.0;
        u = 0;
        /* 999.96 k rounds to 1000.0: move up one unit instead */
        while(u < 5 && floor(x * 10.0 + 0.5) / 10.0 >= 1000.0){
                x /= 1000.0;
                u++;
        }
        fmt_dec(d, sizeof(d), 1, x);
        snprintf(buf, n, "%s %s%s", d, prefix[u], unit);
}

static void fmt_si(char *buf, size_t n, uint64_t v)
{
        fmt_si_unit(buf, n, v, "");
}

/* percent with one decimal and a percent sign, or "n/a" for a value that does not apply */
static void fmt_pct_tile(char *buf, size_t n, double v)
{
        char d[32];
        if(COLLECT_IS_NA(v)){
                snprintf(buf, n, "n/a");
        }else{
                fmt_dec(d, sizeof(d), 1, v);
                snprintf(buf, n, "%s%%", d);
        }
}

static int is_paired(const struct stat_collection *s)
{
        return s->n_read2 > 0;
}

static uint64_t n_mapped(const struct stat_collection *s)
{
        uint64_t n = 0;
        for(int r = 0; r < COLLECT_N_MATE; r++){
                for(int m = 0; m < COLLECT_MAPQ_UNMAPPED; m++){
                        n += s->basic_nums[r][m];
                }
        }
        return n;
}

/* Format as detected by htslib ("BAM", "FASTQ, gzip", ...), taken from the reader's own detection so
 * that standard input is covered too. Empty when there is none or htslib does not recognise it:
 * the header line then leaves the format out. */
static void input_format_name(const htsFormat *f, char *buf, size_t n)
{
        const char *name = NULL;
        int gz;

        buf[0] = 0;
        if(!f){
                return;
        }
        switch(f->format){
        case sam:
                name = "SAM";
                break;
        case bam:
                name = "BAM";
                break;
        case cram:
                name = "CRAM";
                break;
        case fastq_format:
                name = "FASTQ";
                break;
        case fasta_format:
                name = "FASTA";
                break;
        default:
                break;
        }
        gz = (f->compression == gzip || f->compression == bgzf);
        if(name){
                if(gz && f->format != bam && f->format != cram){
                        snprintf(buf, n, "%s, gzip", name);
                }else{
                        snprintf(buf, n, "%s", name);
                }
        }
}

/* ------------------------------------------------------------------ */
/* the page                                                             */
/* ------------------------------------------------------------------ */

int stat_report(struct stat_collection *s, struct samstat_param *p, int id, const htsFormat *fmt)
{
        char *n = NULL;
        tld_strbuf *out = NULL;
        char format[32];

        RUN(tld_get_filename((char *) samstat_output_stem(p, id), &n));

        RUN(stat_collection_finalise(s));
        input_format_name(fmt, format, sizeof(format));

        RUN(tld_strbuf_alloc(&out, 1 << 16));

        RUN(page_head(out, n));
        RUN(app(out, "<body>\n<div class=\"page\">\n"));
        RUN(page_header(out, n, format, s));
        RUN(app(out, "<main>\n"));
        if(s->n_total == 0){
                /* an empty input is valid: say so and leave out the tiles, table and plots */
                RUN(app(out, "<p class=\"banner\">"));
                RUN(esc(out, report_banner_no_reads));
                RUN(app(out, "</p>\n"));
        }else{
                RUN(page_strip(out, s));
                RUN(page_notes(out, s));
                if(s->is_aligned){
                        RUN(page_mapq_table(out, s));
                }
                for(int plot_id = 1; plot_id <= 6; plot_id++){
                        /* plots 5 and 6 exist for aligned input only; nothing is written otherwise */
                        if(plot_id >= 5 && !s->is_aligned){
                                continue;
                        }
                        RUN(page_plot(out, s, plot_id));
                }
        }
        RUN(app(out, "</main>\n"));
        RUN(page_footer(out, n, s));
        RUN(app(out, "</div>\n</body>\n</html>\n"));

        RUN(write_buffer(out, p, id));
        if(p->json){
                RUN(json_report_write(s, p, id));
        }
        if(p->multiqc){
                RUN(mqc_report_write(s, p, id));
        }

        tld_strbuf_free(out);
        MFREE(n);
        return OK;
ERROR:
        if(out){
                tld_strbuf_free(out);
        }
        if(n){
                MFREE(n);
        }
        return FAIL;
}

int write_buffer(tld_strbuf *out_buffer, struct samstat_param *p, int id)
{
        char *buffer = NULL;
        char *filename = NULL;
        struct outfile *o = NULL;
        const char *str = TLD_STR(out_buffer);
        const char *stem = samstat_output_stem(p, id);
        size_t len = strlen(str);

        if(p->outdir){
                RUN(tld_get_filename((char *) stem, &filename));
                RUN(outfile_path(&buffer, "%s/%s.samstat.html", p->outdir, filename));
        }else{
                RUN(outfile_path(&buffer, "%s.samstat.html", stem));
        }
        if(p->verbose){
                LOG_MSG("Writing to: %s", buffer);
        }
        /* written to a temporary file and renamed: a failed run leaves no
           truncated report and never writes through a symlink */
        RUN(outfile_open(&o, buffer));
        ASSERT(fwrite(str, 1, len, o->f) == len, "Failed to write output file: %s", buffer);
        RUN(outfile_commit(&o));
        MFREE(buffer);
        if(filename){
                MFREE(filename);
        }
        return OK;
ERROR:
        outfile_abort(&o);
        if(buffer){
                MFREE(buffer);
        }
        if(filename){
                MFREE(filename);
        }
        return FAIL;
}

/* <head> with the whole style sheet. Light theme only, system fonts, nothing external. */
static const char page_css_1[] =
        ":root{--ink:#1a1a1a;--muted:#555555;--line:#c9cdd3;--wash:#f4f5f7;--accent:#0b4f8a;--accent-bg:#eaf1f8}\n"
        "*{box-sizing:border-box}\n"
        "html{-webkit-text-size-adjust:100%;text-size-adjust:100%}\n"
        "body{margin:0;background:#ffffff;color:var(--ink);"
        "font:16px/1.55 system-ui,-apple-system,\"Segoe UI\",Roboto,\"Helvetica Neue\",Arial,sans-serif}\n"
        ".page{max-width:52rem;margin:0 auto;padding:28px 16px 48px}\n"
        "h1,h2{line-height:1.25;overflow-wrap:anywhere}\n"
        "h1{font-size:1.75rem;margin:0 0 .25rem}\n"
        "h2{font-size:1.3rem;margin:0 0 .35rem}\n"
        "p{margin:0 0 .75rem}\n"
        ".eyebrow{margin:0 0 .25rem;font-size:.875rem;font-weight:600;letter-spacing:.04em;"
        "text-transform:uppercase;color:var(--accent)}\n"
        ".meta{margin:0;color:var(--muted)}\n"
        ".partial{font-weight:600;color:var(--ink)}\n"
        "main>section,main>.strip,main>.notes{margin-top:2.25rem}\n"
        "main>.strip{margin-top:1.75rem}\n"
        /* summary strip */
        ".strip{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin-bottom:0}\n"
        "@media (max-width:42rem){.strip{grid-template-columns:repeat(2,minmax(0,1fr))}}\n"
        "@media (max-width:24rem){.strip{grid-template-columns:minmax(0,1fr)}}\n"
        ".strip>div{margin:0;padding:12px 14px;background:var(--wash);border:1px solid var(--line);border-radius:8px}\n"
        ".strip dt{margin:0;font-size:.9rem;color:var(--muted)}\n"
        ".strip dd{margin:.1rem 0 0;font-size:clamp(1.25rem,5vw,1.75rem);line-height:1.2;font-weight:650;"
        "font-variant-numeric:tabular-nums;overflow-wrap:anywhere}\n"
        ".strip .sub{display:block;font-size:.85rem;font-weight:400;color:var(--muted)}\n"
        /* notes and the banner of an empty file */
        ".notes{list-style:none;padding:0}\n"
        ".banner{margin:1.75rem 0 0;padding:.9rem 1.1rem;background:var(--wash);border-left:5px solid var(--muted);"
        "border-radius:0 6px 6px 0;font-size:1.15rem;font-weight:600}\n"
        ".notes li{margin:0 0 8px;padding:.6rem .9rem;background:var(--wash);border-left:5px solid var(--muted);"
        "border-radius:0 6px 6px 0}\n"
        /* table */
        ".caption{max-width:48rem;margin-bottom:1rem}\n"
        ".table-wrap{overflow-x:auto}\n"
        "table{border-collapse:collapse;min-width:20rem}\n"
        "th,td{padding:.55rem .9rem;border-bottom:1px solid var(--line);text-align:left;vertical-align:middle}\n"
        "thead th{border-bottom:2px solid var(--ink);font-size:.9rem}\n"
        "td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}\n"
        "td svg{display:block}\n"
        ;

/* second half of the style sheet (a string literal is limited to 4095 characters) */
static const char page_css_2[] =
        /* plots */
        ".plot h2 .num{display:inline-block;min-width:1.6rem;color:var(--accent)}\n"
        /* a single figure fills the column; on a narrow screen it keeps 30rem and scrolls sideways so that
           its text stays readable. A pair shares the column side by side, and stacks on a narrow screen. */
        ".figs{margin:0 0 .5rem;overflow-x:auto}\n"
        ".figs svg{display:block;width:100%;margin:0 0 12px 0}\n"
        ".figs:not(.paired)>svg{min-width:30rem}\n"
        ".figs .plot-key svg{width:auto}\n"
        ".figs.paired{max-width:none;display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:0 16px}\n"
        ".figs.paired>p{grid-column:1/-1}\n"
        ".figs.paired svg{max-width:32rem}\n"
        "@media (max-width:40rem){.figs.paired{grid-template-columns:minmax(0,1fr)}}\n"
        ".good{max-width:48rem;margin:0 0 .6rem}\n"
        ".plot-na{display:inline-block;margin:0 0 .75rem;padding:.6rem .9rem;border:2px dashed var(--muted);border-radius:6px}\n"
        ".plot-note{margin:0 0 .75rem;padding:.5rem .9rem;background:var(--wash);border-left:5px solid var(--muted);"
        "border-radius:0 6px 6px 0;max-width:48rem}\n"
        ".plot-key{margin:0 0 .5rem}\n"
        ".plot-key svg{margin:0}\n"
        "details{margin:.25rem 0 0;border:1px solid var(--line);border-radius:6px;max-width:48rem}\n"
        "summary{padding:.5rem .9rem;font-weight:600;color:var(--accent);cursor:pointer}\n"
        "summary:focus-visible{outline:3px solid var(--accent);outline-offset:2px}\n"
        "details[open] summary{border-bottom:1px solid var(--line)}\n"
        ".howto{padding:.75rem .9rem 0}\n"
        ".howto p{margin-bottom:.6rem}\n"
        /* footer */
        "footer{margin-top:3rem;padding-top:1rem;border-top:2px solid var(--line);font-size:.9rem;color:var(--muted)}\n"
        "footer p{margin-bottom:.5rem}\n"
        /* print: A4 portrait, every explanation open, figures stay whole */
        "@page{size:A4 portrait;margin:14mm}\n"
        "@media print{\n"
        "body{font-size:11pt;-webkit-print-color-adjust:exact;print-color-adjust:exact}\n"
        ".page{max-width:none;padding:0}\n"
        "main>section,.strip>div,.notes li,details,table,tr,.figs,.plot-na,.plot-note,.good{break-inside:avoid}\n"
        "h1,h2,summary{break-after:avoid}\n"
        "section.plot{break-inside:avoid-page}\n"
        "details>*:not(summary){display:block!important}\n"
        "details::details-content{content-visibility:visible!important;display:block!important}\n"
        "details>summary{list-style:none;border-bottom:1px solid var(--line)}\n"
        "details>summary::-webkit-details-marker{display:none}\n"
        "summary{color:var(--ink)}\n"
        "}\n";

static int page_head(tld_strbuf *o, const char *name)
{
        RUN(app(o, "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n"));
        RUN(app(o, "<meta charset=\"utf-8\">\n"));
        RUN(app(o, "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"));
        RUN(app(o, "<meta name=\"color-scheme\" content=\"light\">\n"));
        RUN(app(o, "<title>"));
        RUN(esc(o, name));
        RUN(app(o, " - SAMStat report</title>\n"));
        RUN(app(o, "<style>\n"));
        RUN(app(o, page_css_1));
        RUN(app(o, page_css_2));
        RUN(app(o, "</style>\n</head>\n"));
        return OK;
ERROR:
        return FAIL;
}

/* 1.1: <file> - <format> - <layout> - SAMStat <version> [- partial] */
static int page_header(tld_strbuf *o, const char *name, const char *format, struct stat_collection *s)
{
        char buf[MAX_NUM];
        char len[MAX_NUM];

        RUN(app(o, "<header>\n<p class=\"eyebrow\">SAMStat report</p>\n<h1>"));
        RUN(esc(o, name));
        RUN(app(o, "</h1>\n<p class=\"meta\">"));
        if(format[0]){
                RUN(esc(o, format));
                RUN(app(o, " &middot; "));
        }
        if(s->n_total && s->total_bases){
                fmt_len(len, sizeof(len), s->median_best, s->median_exact);
                RUN(app(o, s->layout == COLLECT_LAYOUT_LONG ? "long" : "short"));
                RUN(app(o, " reads (median "));
                RUN(esc(o, len));
                RUN(app(o, ") &middot; "));
        }
        RUN(app(o, "SAMStat "));
        RUN(esc(o, PACKAGE_VERSION));
        if(s->is_partial_report){
                fmt_count(buf, sizeof(buf), s->n_total);
                RUN(app(o, "<span class=\"partial\">"));
                RUN(app_template(o, report_header_partial, buf, NULL));
                RUN(app(o, "</span>"));
        }
        RUN(app(o, "</p>\n</header>\n"));
        return OK;
ERROR:
        return FAIL;
}

/* 1.2: one card per value, same order every time */
static int strip_item(tld_strbuf *o, const char *label, const char *value, const char *sub)
{
        RUN(app(o, "<div><dt>"));
        RUN(esc(o, label));
        RUN(app(o, "</dt><dd>"));
        RUN(esc(o, value));
        if(sub && sub[0]){
                RUN(app(o, "<span class=\"sub\">"));
                RUN(esc(o, sub));
                RUN(app(o, "</span>"));
        }
        RUN(app(o, "</dd></div>\n"));
        return OK;
ERROR:
        return FAIL;
}

static int page_strip(tld_strbuf *o, struct stat_collection *s)
{
        char v[MAX_NUM];
        char sub[MAX_NUM];
        char tmp[MAX_NUM];
        const int has_seq = s->total_bases > 0;

        RUN(app(o, "<dl class=\"strip\" aria-label=\"Summary\">\n"));

        fmt_si(v, sizeof(v), s->n_total);
        sub[0] = 0;
        if(s->n_total >= 10000){
                fmt_count(sub, sizeof(sub), s->n_total);
        }
        RUN(strip_item(o, is_paired(s) ? "Reads (R1+R2)" : "Reads", v, sub));

        /* no sequence (SEQ is * in every read): length values do not exist, they are not 0 */
        sub[0] = 0;
        if(has_seq){
                fmt_si_unit(v, sizeof(v), s->total_bases, "b");
                if(s->total_bases >= 10000){
                        fmt_count(sub, sizeof(sub), s->total_bases);
                }
        }else{
                snprintf(v, sizeof(v), "n/a");
        }
        RUN(strip_item(o, "Bases", v, sub));

        if(has_seq){
                fmt_len(v, sizeof(v), s->n50_best, s->n50_exact);
        }else{
                snprintf(v, sizeof(v), "n/a");
        }
        RUN(strip_item(o, "N50", v, has_seq ? report_tile_n50_sub : NULL));

        /* mean Q is left out when no base has a quality (FASTA) */
        if(s->qual_bases_all > 0){
                fmt_dec(v, sizeof(v), 1, s->mean_q);
                RUN(strip_item(o, "Mean Q", v, report_tile_meanq_sub));
        }

        if(s->is_aligned){
                fmt_pct_tile(v, sizeof(v), s->pct_mapped);
                RUN(strip_item(o, "% mapped", v, NULL));
                fmt_pct_tile(v, sizeof(v), s->pct_mapq_ge30);
                RUN(strip_item(o, "% MAPQ >= 30", v, COLLECT_IS_NA(s->pct_mapq_ge30) ? NULL : report_tile_mapq_sub));
                fmt_pct_tile(v, sizeof(v), s->pct_clipped);
                RUN(strip_item(o, "% bases clipped", v, COLLECT_IS_NA(s->pct_clipped) ? NULL : report_tile_clipped_sub));
                if(COLLECT_IS_NA(s->median_identity)){
                        snprintf(v, sizeof(v), "n/a");
                }else{
                        fmt_dec(tmp, sizeof(tmp), 1, s->median_identity);
                        snprintf(v, sizeof(v), "%s%%", tmp);
                }
                RUN(strip_item(o, "Median identity", v, COLLECT_IS_NA(s->median_identity) ? NULL : report_tile_identity_sub));
        }
        RUN(app(o, "</dl>\n"));
        return OK;
ERROR:
        return FAIL;
}

/* 1.3: facts, never verdicts. Only the lines whose condition holds. */
static int note_item(tld_strbuf *o, const char *tmpl, const char *a1, const char *a2)
{
        RUN(app(o, "<li>"));
        RUN(app_template(o, tmpl, a1, a2));
        RUN(app(o, "</li>\n"));
        return OK;
ERROR:
        return FAIL;
}

static int page_notes(tld_strbuf *o, struct stat_collection *s)
{
        char a[MAX_NUM];
        char b[MAX_NUM];
        uint64_t mapped = n_mapped(s);
        int partial = s->is_partial_report;
        int no_md = s->reads_without_md > 0;
        int ident_na = s->identity_na > 0;
        int beyond = (s->layout == COLLECT_LAYOUT_SHORT && s->bases_beyond_position_axis > 0);
        int no_seq = (s->n_total > 0 && s->total_bases == 0);
        int hidden = (!no_seq && s->n_total > 0 && s->n_total < REPORT_MIN_POS_BASES);   /* plots 3-5 are not drawn */
        int few = (!no_seq && !hidden && s->n_total > 0 && s->n_total < REPORT_FEW_READS);
        uint64_t in_axis = 0;
        double pct = 0.0;

        if(!(partial || no_md || ident_na || beyond || no_seq || few || hidden)){
                return OK;
        }
        RUN(app(o, "<ul class=\"notes\" aria-label=\"Notes\">\n"));
        if(partial){
                fmt_count(a, sizeof(a), s->n_total);
                RUN(note_item(o, report_note_partial, a, NULL));
        }
        if(no_seq){
                RUN(note_item(o, report_note_no_sequence, NULL, NULL));
        }
        if(hidden){
                fmt_count(a, sizeof(a), s->n_total);
                RUN(note_item(o, report_note_hidden, a, s->n_total == 1 ? "" : "s"));
        }
        if(few){
                fmt_count(a, sizeof(a), s->n_total);
                RUN(note_item(o, report_note_few_reads, a, REPORT_MIN_POS_BASES_STR));
        }
        if(no_md){
                fmt_count(a, sizeof(a), s->reads_without_md);
                fmt_count(b, sizeof(b), mapped);
                RUN(note_item(o, report_note_no_md, a, b));
        }
        if(ident_na){
                fmt_count(a, sizeof(a), s->identity_na);
                fmt_count(b, sizeof(b), mapped);
                RUN(note_item(o, report_note_identity_na, a, b));
        }
        if(beyond){
                /* share of the bases of the position plots that lie past the last bin */
                for(int r = 0; r < COLLECT_N_MATE; r++){
                        for(int m = 0; m < COLLECT_N_MAPQ; m++){
                                for(int i = 0; i < COLLECT_ABS_BINS; i++){
                                        const struct collect_pos_bin *pb = &s->acc->pos_abs[r][m][i];
                                        for(int k = 0; k < COLLECT_N_BASE; k++){
                                                in_axis += pb->comp[k];
                                        }
                                }
                        }
                }
                pct = 100.0 * (double) s->bases_beyond_position_axis
                        / (double) (s->bases_beyond_position_axis + in_axis);
                if(pct > 0.0 && pct < 0.1){
                        /* template arguments are written as they are, so the "<" is escaped here */
                        snprintf(a, sizeof(a), "&lt;0.1");
                }else{
                        fmt_dec(a, sizeof(a), 1, pct);
                }
                RUN(note_item(o, report_note_beyond_axis, a, NULL));
        }
        RUN(app(o, "</ul>\n"));
        return OK;
ERROR:
        return FAIL;
}

/* 1.4: the table that is also the key of every MAPQ line on the page */
static int page_mapq_table(tld_strbuf *o, struct stat_collection *s)
{
        char v[MAX_NUM];
        int paired = is_paired(s);

        RUN(app(o, "<section id=\"mapq\" aria-labelledby=\"mapq-title\">\n<h2 id=\"mapq-title\">"));
        RUN(esc(o, report_mapq_title));
        RUN(app(o, "</h2>\n<p class=\"caption\">"));
        RUN(esc(o, report_mapq_caption));
        RUN(app(o, "</p>\n<div class=\"table-wrap\">\n<table>\n<thead><tr>"));
        RUN(app(o, "<th scope=\"col\">Line in plots</th><th scope=\"col\">Class</th>"));
        RUN(app(o, "<th scope=\"col\" class=\"num\">Reads</th><th scope=\"col\" class=\"num\">% of reads</th>"));
        if(paired){
                RUN(app(o, "<th scope=\"col\" class=\"num\">Read 1</th><th scope=\"col\" class=\"num\">Read 2</th>"));
        }
        RUN(app(o, "</tr></thead>\n<tbody>\n"));
        for(int m = 0; m < COLLECT_N_MAPQ; m++){
                uint64_t reads = s->basic_nums[0][m] + s->basic_nums[1][m];
                RUN(app(o, "<tr><td>"));
                RUN(report_mapq_swatch(o, m));
                RUN(app(o, "</td><th scope=\"row\">"));
                RUN(esc(o, report_mapq_help[m].label));
                RUN(app(o, "</th><td class=\"num\">"));
                fmt_count(v, sizeof(v), reads);
                RUN(app(o, v));
                RUN(app(o, "</td><td class=\"num\">"));
                if(s->n_total){
                        fmt_dec(v, sizeof(v), 1, 100.0 * (double) reads / (double) s->n_total);
                }else{
                        snprintf(v, sizeof(v), "n/a");
                }
                RUN(app(o, v));
                RUN(app(o, "</td>"));
                if(paired){
                        for(int r = 0; r < COLLECT_N_MATE; r++){
                                RUN(app(o, "<td class=\"num\">"));
                                fmt_count(v, sizeof(v), s->basic_nums[r][m]);
                                RUN(app(o, v));
                                RUN(app(o, "</td>"));
                        }
                }
                RUN(app(o, "</tr>\n"));
        }
        RUN(app(o, "</tbody>\n</table>\n</div>\n</section>\n"));
        return OK;
ERROR:
        return FAIL;
}

/* One plot: question, one-sentence caption, the figure(s), the "What good looks like:" line, collapsed explanation. */
static int page_plot(tld_strbuf *o, struct stat_collection *s, int plot_id)
{
        const struct report_plot_help *h = &report_plot_help[plot_id - 1];
        char num[16];
        tld_strbuf *fig = NULL;
        int has_figure;

        snprintf(num, sizeof(num), "%d", plot_id);
        RUN(app(o, "<section class=\"plot\" id=\"plot-"));
        RUN(app(o, num));
        RUN(app(o, "\" aria-labelledby=\"plot-"));
        RUN(app(o, num));
        RUN(app(o, "-title\">\n<h2 id=\"plot-"));
        RUN(app(o, num));
        RUN(app(o, "-title\"><span class=\"num\">"));
        RUN(app(o, num));
        RUN(app(o, ".</span> "));
        RUN(esc(o, h->title));
        RUN(app(o, "</h2>\n"));
        /* rendered first: a plot can be replaced by one sentence (one read length, too few reads), and then
         * the caption and the explanation would describe a figure that is not there */
        RUN(tld_strbuf_alloc(&fig, 1 << 14));
        RUN(report_plot_render(fig, s, plot_id));
        has_figure = (strstr(TLD_STR(fig), "<svg") != NULL);
        if(has_figure){
                RUN(app(o, "<p class=\"caption\">"));
                RUN(esc(o, report_plot_caption(plot_id, s->is_aligned)));
                RUN(app(o, "</p>\n"));
        }
        RUN(app(o, is_paired(s) && has_figure ? "<div class=\"figs paired\">\n" : "<div class=\"figs\">\n"));
        RUN(app(o, TLD_STR(fig)));
        tld_strbuf_free(fig);
        fig = NULL;
        RUN(app(o, "</div>\n"));

        if(has_figure){
                RUN(app(o, "<p class=\"good\"><strong>"));
                RUN(esc(o, report_good_label));
                RUN(app(o, "</strong> "));
                RUN(esc(o, report_plot_good(plot_id, s->is_aligned)));
                RUN(app(o, "</p>\n"));
                /* progressive disclosure: the counts above are the message; the shape view is opened on purpose */
                RUN(tld_strbuf_alloc(&fig, 1 << 14));
                RUN(report_plot_render_shape(fig, s, plot_id));
                if(strstr(TLD_STR(fig), "<svg") != NULL){
                        RUN(app(o, "<details class=\"shape\">\n<summary>"));
                        RUN(esc(o, report_shape_label));
                        RUN(app(o, "</summary>\n<div class=\"howto\">\n<p>"));
                        RUN(esc(o, report_shape_caption));
                        RUN(app(o, "</p>\n"));
                        RUN(app(o, is_paired(s) ? "<div class=\"figs paired\">\n" : "<div class=\"figs\">\n"));
                        RUN(app(o, TLD_STR(fig)));
                        RUN(app(o, "</div>\n</div>\n</details>\n"));
                }
                tld_strbuf_free(fig);
                fig = NULL;
                RUN(app(o, "<details>\n<summary>"));
                RUN(esc(o, report_how_to_read_label));
                RUN(app(o, " "));
                RUN(app(o, num));
                RUN(app(o, "</summary>\n<div class=\"howto\">\n"));
                RUN(app(o, "<p><strong>"));
                RUN(esc(o, report_look_label));
                RUN(app(o, "</strong> "));
                RUN(esc(o, report_plot_look(plot_id, s->is_aligned)));
                RUN(app(o, "</p>\n<p><strong>"));
                RUN(esc(o, report_benign_label));
                RUN(app(o, "</strong> "));
                RUN(esc(o, h->benign));
                RUN(app(o, "</p>\n"));
                if(h->note){
                        RUN(app(o, "<p><strong>"));
                        RUN(esc(o, report_note_label));
                        RUN(app(o, "</strong> "));
                        RUN(esc(o, h->note));
                        RUN(app(o, "</p>\n"));
                }
                RUN(app(o, "</div>\n</details>\n"));
        }
        RUN(app(o, "</section>\n"));
        return OK;
ERROR:
        tld_strbuf_free(fig);
        return FAIL;
}

/* 1.6 / F10: version, input basename, citation, raw secondary and supplementary counts.
 * No command line (it would leak absolute paths and differ with -t). Plain text only: the
 * report makes no request and links nowhere. */
static int page_footer(tld_strbuf *o, const char *name, struct stat_collection *s)
{
        char v[MAX_NUM];

        RUN(app(o, "<footer>\n<p>Report of <strong>"));
        RUN(esc(o, name));
        RUN(app(o, "</strong> made by "));
        RUN(esc(o, PACKAGE_NAME));
        RUN(app(o, " "));
        RUN(esc(o, PACKAGE_VERSION));
        RUN(app(o, ".</p>\n"));
        if(s->n_secondary || s->n_supplementary){
                RUN(app(o, "<p>"));
                if(s->n_secondary){
                        fmt_count(v, sizeof(v), s->n_secondary);
                        RUN(app(o, v));
                        RUN(app(o, " secondary alignments (not counted)"));
                }
                if(s->n_secondary && s->n_supplementary){
                        RUN(app(o, "; "));
                }
                if(s->n_supplementary){
                        fmt_count(v, sizeof(v), s->n_supplementary);
                        RUN(app(o, v));
                        RUN(app(o, " supplementary alignments (not counted)"));
                }
                RUN(app(o, ".</p>\n"));
        }
        RUN(app(o, "<p>Please cite: Lassmann T., &ldquo;SAMStat 2: quality control for next generation sequencing data&rdquo;, "
                "<em>Bioinformatics</em> 39(1), btad019, 2023, doi:10.1093/bioinformatics/btad019.</p>\n"));
        RUN(app(o, "<p>Original method: Lassmann T., Hayashizaki Y., Daub C.O., &ldquo;SAMStat: monitoring biases in next "
                "generation sequencing data&rdquo;, <em>Bioinformatics</em> 27(1), 130-131, 2011, "
                "doi:10.1093/bioinformatics/btq614.</p>\n"));
        RUN(app(o, "</footer>\n"));
        return OK;
ERROR:
        return FAIL;
}
