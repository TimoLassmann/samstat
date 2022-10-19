#ifndef PLOT_INTERVAL_H
#define PLOT_INTERVAL_H

#ifdef PLOT_INTERVAL_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int get_plot_interval(int len,int bin_start,int target_n_bin, int *interval);

#undef PLOT_INTERVAL_IMPORT
#undef EXTERN


#endif
