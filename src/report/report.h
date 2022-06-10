#ifndef REPORT_H
#define REPORT_H


#include "../metrics/metrics.h"
#include "../param/param.h"

#ifdef REPORT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int create_report(struct metrics*m, struct samstat_param* p);



#undef REPORT_IMPORT
#undef EXTERN


#endif
