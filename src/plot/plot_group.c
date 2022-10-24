#include "tld.h"

#include "plot.h"


#define PLOT_GROUP_IMPORT
#include "plot_group.h"



static int set_series_plot(struct plot_data *d,struct plot_group *g);
static int debug_print_group_data(struct plot_group *g);
static int sanity_check_group_data(struct plot_group *g);
static int set_visibility(struct plot_group *g);
static int check_for_group_data(struct plot_data *d, struct plot_group* g);

int process_into_groups(struct plot_data *d, struct plot_group** groups)
{
        struct plot_group* g = NULL;
        int n_group = 0;
        /* d->n_group = 0; */

        if(!d->group_size || d->L < d->group_size){
                RUN(plot_groups_alloc(&g, 1));
                g->l[0]->start = 0;
                g->l[0]->stop = d->L;
                g->l[0]->plot = 1;
                g->l[0]->viz = 1;
        }else if(d->group_size){
                int start = 0;
                int g_id = 0;
                while(start < d->L){

                        start += d->group_size;
                        n_group++;
                }
                RUN(plot_groups_alloc(&g, n_group));
                start = 0;
                while(start < d->L){
                        g->l[g_id]->start = start;
                        g->l[g_id]->stop = start + d->group_size;
                        g->l[g_id]->viz = 0;
                        if(g_id == 0){
                                g->l[g_id]->viz = 1;
                        }
                        start += d->group_size;
                        g_id++;
                }
        }else{
                ERROR_MSG("This should never happen");
        }

        /* Make sure there is data  */
        RUN(check_for_group_data(d, g));

        /* Sanity checks  */
        sanity_check_group_data(g);
        /* Set series to plot or not */
        set_series_plot(d, g);
        /* Make sure visibility is set correctly */
        set_visibility(g);
        /* Sanity checks  */
        sanity_check_group_data(g);

        /* Debug print */
        /* debug_print_group_data(g); */

        *groups = g;
        return OK;
ERROR:
        plot_groups_free(g);
        return FAIL;
}

int set_series_plot(struct plot_data *d,struct plot_group *g)
{
        ASSERT(g != NULL,"No group data");

        ASSERT(g->n_groups > 0, "No groups !");
        struct plot_group_itm* itm = NULL;
        int plot_empty  = d->plot_empty_series;


        for(int i = 0; i < g->n_groups;i++){
                itm = g->l[i];
                int s = itm->start;
                int e = itm->stop;
                galloc(&itm->plot_series, d->L);
                for(int j = 0; j < d->L;j++){
                        itm->plot_series[j] = 0;
                }
                for(int j = s; j < e;j++){
                        itm->plot_series[j] = 1;
                }
        }
        if(plot_empty == 0){
                uint64_t** l_d = NULL;
                int l_len;
                l_d = d->data;
                l_len = d->len;
                if(d->clu_len){
                        l_d = d->clu_data;
                        l_len = d->clu_len;
                }

                for(int i = 0; i < g->n_groups;i++){
                        itm = g->l[i];
                        int s = itm->start;
                        int e = itm->stop;

                        for(int j = s; j < e;j++){
                                uint64_t sum = 0;
                                for(int c = 0; c < l_len;c++){
                                        if(l_d[j][c]){
                                                sum = 1;
                                                break;
                                        }
                                }
                                if(!sum){
                                        itm->plot_series[j] = 0;
                                }
                        }
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int set_visibility(struct plot_group *g)
{

        ASSERT(g != NULL,"No group data");

        ASSERT(g->n_groups > 0, "No groups !");
        int vis = 1;
        for(int i = 0; i < g->n_groups;i++){
                g->l[i]->viz = 0;
                if(g->l[i]->plot){
                        g->l[i]->viz = vis;
                        vis = 0;
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int sanity_check_group_data(struct plot_group *g)
{
        ASSERT(g != NULL,"No group data");

        ASSERT(g->n_groups > 0, "No groups !");

        for(int i = 0; i < g->n_groups;i++){
                if(g->l[i]->start >= g->l[i]->stop){
                        ERROR_MSG("in group %d: start is >= stop: %d %d", g->l[i]->start, g->l[i]->stop);
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int debug_print_group_data(struct plot_group *g)
{
        ASSERT(g != NULL,"No group data");

        ASSERT(g->n_groups > 0, "No groups !");


        LOG_MSG("Will plot: %d", g->plot);
        for(int i = 0; i < g->n_groups;i++){
                LOG_MSG("Group%4d: %d %d %d %d",i,g->l[i]->start,g->l[i]->stop, g->l[i]->plot, g->l[i]->viz);
        }
        return OK;
ERROR:
        return FAIL;
}


int check_for_group_data(struct plot_data *d, struct plot_group* g)
{
        uint64_t** l_d = NULL;
        int l_len;
        l_d = d->data;
        l_len = d->len;
        if(d->clu_len){
                l_d = d->clu_data;
                l_len = d->clu_len;
        }

        for(int i = 0; i < g->n_groups;i++){
                int s = g->l[i]->start;
                int e = g->l[i]->stop;
                uint64_t sum = 0;
                for(int j = s; j < e;j++){
                        for(int c = 0; c < l_len;c++){
                                if(l_d[j][c]){
                                        sum = 1;
                                }
                        }
                }
                g->l[i]->plot = 1;
                if(!sum){
                        g->l[i]->plot = 0;
                        g->l[i]->viz = 0;
                }
        }
        g->plot = 0;
        for(int i = 0; i < g->n_groups;i++){
                if(g->l[i]->plot){
                        g->plot = 1;
                }
        }
        return OK;
}

int plot_groups_alloc(struct plot_group **plot_group, int n)
{
        struct plot_group* g = NULL;
        ASSERT(n > 0,"There must be at least one plot group");

        MMALLOC(g, sizeof(struct plot_group));

        g->n_groups = n;
        g->plot = 1;
        MMALLOC(g->l, sizeof(struct plot_group_itm*) * g->n_groups);
        for(int i = 0; i < g->n_groups;i++){
                g->l[i] = NULL;
                MMALLOC(g->l[i], sizeof(struct plot_group_itm));
                g->l[i]->start = 0;
                g->l[i]->stop = 0;
                g->l[i]->plot = 1;
                g->l[i]->viz = 1;
                g->l[i]->plot_series = NULL;
        }

        *plot_group = g;

        return OK;
ERROR:
        plot_groups_free(g);
        return FAIL;
}

void plot_groups_free(struct plot_group *g)
{
        if(g){
                if(g->n_groups){
                        for(int i = 0; i < g->n_groups;i++){
                                if(g->l[i]->plot_series){
                                        gfree(g->l[i]->plot_series);
                                }
                                MFREE(g->l[i]);
                        }
                        MFREE(g->l);
                }
                MFREE(g);
        }
}
