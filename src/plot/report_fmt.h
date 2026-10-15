#ifndef REPORT_FMT_H
#define REPORT_FMT_H

/*
 * Number formatting shared by the report page (report/stat_report.c) and the plots
 * (report_plots.c), so that a count or a length reads the same in the summary strip and
 * in a plot. Header only (static inline): no extra source file to build, and a file that
 * uses only some of the helpers gets no unused-function warning.
 *
 * Every helper writes a NUL-terminated string of at most n bytes into buf and does not
 * depend on the locale.
 */

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* 1234567 -> "1,234,567" */
static inline void fmt_count(char *buf, size_t n, uint64_t v)
{
        char raw[32];
        char out[48];
        int len = snprintf(raw, sizeof(raw), "%" PRIu64, v);
        int k = 0;
        for(int i = 0; i < len; i++){
                if(i > 0 && (len - i) % 3 == 0){
                        out[k++] = ',';
                }
                out[k++] = raw[i];
        }
        out[k] = 0;
        snprintf(buf, n, "%s", out);
}

/* "%.<decimals>f" that does not depend on the locale */
static inline void fmt_dec(char *buf, size_t n, int decimals, double v)
{
        snprintf(buf, n, "%.*f", decimals, v);
        for(char *c = buf; *c; c++){
                if(*c == ','){
                        *c = '.';
                }
        }
}

/* Lengths: bp below 10 kb, kb up to 1 Mb, Mb above (the unit of the Bases tile). Above 1 kb a value
 * that is a bin edge (within 1%) is marked with "~"; exact is set when the value is the true length. */
static inline void fmt_len(char *buf, size_t n, uint64_t v, int exact)
{
        char t[48];
        const char *tilde = (v < 1000 || exact) ? "" : "~";
        if(v < 1000){
                fmt_count(t, sizeof(t), v);
                snprintf(buf, n, "%s bp", t);
        }else if(v < 10000){
                fmt_count(t, sizeof(t), v);
                snprintf(buf, n, "%s%s bp", tilde, t);
        }else if(v >= 999950){
                char d[32];
                fmt_dec(d, sizeof(d), 1, (double) v / 1000000.0);
                snprintf(buf, n, "%s%s Mb", tilde, d);
        }else{
                char d[32];
                char ip[32];
                char *dot;
                size_t il;
                fmt_dec(d, sizeof(d), 1, (double) v / 1000.0);
                /* thousands separators in the integer part of the kb value */
                dot = strchr(d, '.');
                il = dot ? (size_t) (dot - d) : strlen(d);
                snprintf(ip, sizeof(ip), "%.*s", (int) il, d);
                fmt_count(t, sizeof(t), (uint64_t) strtoull(ip, NULL, 10));
                snprintf(buf, n, "%s%s%s kb", tilde, t, dot ? dot : "");
        }
}

#endif
