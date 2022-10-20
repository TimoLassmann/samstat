#ifndef PLOT_INTERVAL_H
#define PLOT_INTERVAL_H

#include <stdint.h>

#ifdef PLOT_INTERVAL_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


struct plot_data;

EXTERN int clu_data(struct plot_data *d);
#undef PLOT_INTERVAL_IMPORT
#undef EXTERN


#endif
