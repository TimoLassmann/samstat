#ifndef PLOT_GROUP_H
#define PLOT_GROUP_H


#include <inttypes.h>
#ifdef PLOT_GROUP_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif
struct plot_data;
struct plot_group_itm {
        uint8_t* plot_series;
        int start;
        int stop;
        int viz;
        int plot;
};


struct plot_group {
        struct plot_group_itm** l;
        int n_groups;
        int plot;
};


EXTERN int process_into_groups(struct plot_data *d, struct plot_group** groups);
EXTERN int plot_groups_alloc(struct plot_group **plot_group, int n);
EXTERN void plot_groups_free(struct plot_group *g);

#undef PLOT_GROUP_IMPORT
#undef EXTERN


#endif
