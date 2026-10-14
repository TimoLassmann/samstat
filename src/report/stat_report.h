#ifndef STAT_REPORT_H
#define STAT_REPORT_H

#ifdef STAT_REPORT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct stat_collection;
struct samstat_param;
struct htsFormat;

EXTERN int stat_report(struct stat_collection* s, struct samstat_param *p, int id, const struct htsFormat *fmt);   /* fmt: as detected by the reader, may be NULL */


#undef STAT_REPORT_IMPORT
#undef EXTERN


#endif
