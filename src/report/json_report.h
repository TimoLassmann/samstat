#ifndef JSON_REPORT_H
#define JSON_REPORT_H

#ifdef JSON_REPORT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct stat_collection;
struct samstat_param;

/* Writes <basename>.samstat.json next to the HTML report.
   Must be called after stat_collection_finalise(). */
EXTERN int json_report_write(struct stat_collection* s, struct samstat_param *p, int id);

#undef JSON_REPORT_IMPORT
#undef EXTERN

#endif
