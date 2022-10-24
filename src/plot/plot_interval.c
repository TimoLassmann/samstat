#include "tld.h"

#include "plot.h"

#define PLOT_INTERVAL_IMPORT
#include "plot_interval.h"

struct interval_clu {
        uint64_t** data;
        double** n;
        char** x_axis_labels;
        int* mapping;
        int L;
        int n_clu;
        int clu_len;
        int org_len;
};

static int get_plot_interval_lin(struct plot_data *pd, struct interval_clu **interval_clu);
static int get_plot_interval_exp(struct plot_data* pd, struct interval_clu** interval_clu);
static int interval_clu_alloc(struct interval_clu **interval_clu, int len);

static int summarise_data(struct plot_data *d, struct interval_clu *ic);

static void interval_clu_free(struct interval_clu *ic);

/*
  Inspired by code from fastqc and specifically:

  https://github.com/s-andrews/FastQC/blob/bfde9d31aa94f29e9c58b6e91574a08e57335a4a/uk/ac/babraham/FastQC/Graphs/BaseGroup.java

*/

int clu_data(struct plot_data *d)
{
        struct interval_clu* ic = NULL;
        ASSERT(d != NULL, "No data!");
        if(d->len <= 250){
                d->x_is_categorical = 0;
                galloc(&d->x_axis_labels, d->len, 32);
                for(int i = 0; i < d->len;i++){
                        snprintf(d->x_axis_labels[i],32,"%d",i+1);
                }
                /* d->clu_data = d->data; */
                /* d->clu_len = d->len; */
        }else if(d->len <= 1000){
                /* needs to be removed ... */
                /* d->bin_start = 50; */
                /* d->target_n_clu = 500; */
                get_plot_interval_lin(d, &ic);
                d->x_is_categorical = 1;
        }else{
                /* d->bin_start = 50; */
                get_plot_interval_exp(d, &ic);
                d->x_is_categorical = 1;
        }
        if(ic){
                RUN(summarise_data(d,ic));
                interval_clu_free(ic);
        }
        return OK;
ERROR:
        if(ic){
                interval_clu_free(ic);
        }
        return FAIL;
}


int summarise_data(struct plot_data *d, struct interval_clu *ic)
{
        /* first need to allocate storage; */
        ic->L = d->L;
        galloc(&ic->data, ic->L, ic->clu_len);
        galloc(&ic->n,  ic->L,ic->clu_len);
        for(int i = 0; i < ic->L;i++){
                for(int j = 0; j < ic->clu_len;j++){
                        ic->data[i][j] = 0;
                        ic->n[i][j] = 0;
                }
        }
        for(int i = 0; i < d->L ;i++){
                for(int j = 0; j < d->len;j++){
                        int idx = ic->mapping[j];
                        /* fprintf(stderr,"%d %d  %d\n",i,j,ic->clu_len); */
                        ic->data[i][idx] += d->data[i][j];
                        ic->n[i][idx] += 1.0;
                }
        }
        /* fprintf(stderr,"Summarise \n"); */
        for(int i = 0; i < d->L ;i++){
                for(int j = 0; j < ic->clu_len;j++){
                        ic->data[i][j] = (uint64_t) (round((double)ic->data[i][j] / ic->n[i][j]));
                        /* fprintf(stderr,"%3ld ",ic->data[i][j] ); */
                }
                /* fprintf(stderr,"\n"); */
        }
        /* fprintf(stderr,"\n"); */
        /* gfree(d->data); */

        d->clu_data = ic->data;
        ic->data = NULL;
        /* d->x_axis_labels = ic->x_axis_labels; */
        d->clu_len = ic->clu_len;
        return OK;
ERROR:
        return FAIL;
}

int get_plot_interval_lin( struct plot_data* pd, struct interval_clu** interval_clu)
{
        struct interval_clu* ic = NULL;
        int target_n_bin = pd->target_n_clu;
        int base[2] = {5,10};
        int ival = 0;
        int igroup = 0;
        int m = 1;
        int run = 1;

        int len = pd->len;
        int bin_start = pd->bin_start;
        int clu_id = 0;
        int old_start = 0;
        interval_clu_alloc(&ic, len);

        /* fprintf(stderr,"LEN: %d %d\n", len,bin_start); */
        while(run){
                for(int i = 0; i < 2;i++){
                        ival = base[i] * m;
                        igroup = bin_start + ((len-bin_start) / ival);
                        if((len-bin_start) % ival != 0){
                                igroup += 1;
                        }
                        /* fprintf(stderr,"group count :  %d (int: %d) all : %d targer n_greoup: %d\n", igroup,ival, bin_start + (igroup - bin_start)*ival ,target_n_bin ); */
                        if(igroup <  target_n_bin){
                                run = 0;
                                break;
                        }
                }
                m *= 10;
                if (m == 10000000) {
                        ERROR_MSG("No interval found");
                }
        }
        /* fprintf(stderr,"val: %d\n", ival); */
        clu_id = 0;

        for(int i = 0; i < len;i++){
                if(i <= bin_start){
                        ic->mapping[i] = clu_id;
                        clu_id++;
                }else{
                        /* clu_id = (i-1) / ival + bin_start -1; */
                        ic->mapping[i] = clu_id;
                        if((i-bin_start) % ival == 0){
                                clu_id++;
                        }
                }
                /* fprintf(stderr,"%d -> %d\n", i, ic->mapping[i]); */
        }

        /* fprintf(stderr,"%4d - %4d\n",old_start,len); */

        ic->clu_len = ic->mapping[len-1] + 1;

        galloc(&pd->x_axis_labels, ic->clu_len, 32);
        old_start = 0;
        clu_id = 0;
        for(int i = 0; i < len;i++){
                if(i <= bin_start){
                        snprintf(pd->x_axis_labels[clu_id],32,"%d",i);
                        old_start = i;
                        clu_id++;
                }else{
                        if((i-bin_start) % ival == 0){
                                /* fprintf(stderr,"%4d - %4d\n",old_start,i); */
                                snprintf(pd->x_axis_labels[clu_id],32,"%d-%d",old_start,i);
                                old_start = i;
                                clu_id++;
                        }
                }
        }
        /* for(int i = 0; i < ic->clu_len;i++){ */
        /*         fprintf(stderr,"%s\n",pd->x_axis_labels[i]); */
        /* } */
        *interval_clu = ic;

        return OK;
ERROR:
        return FAIL;
}


int get_plot_interval_exp(struct plot_data* pd, struct interval_clu** interval_clu)
{

        struct interval_clu* ic = NULL;
        int len = pd->len;
        int bin_start = pd->bin_start;
        int clu_id = 0;
        int old_start = 0;
        int ival = 1;
        interval_clu_alloc(&ic, len);



        /* Vector<BaseGroup> groups = new Vector<BaseGroup>(); */

        clu_id = 0;
        for(int i = 0; i < len;i++){
                if(i <= bin_start){
                        ic->mapping[i] = clu_id;

                        clu_id++;
                }else{

                        if(i < 2000){
                                ival = 10;
                        }else if(i < 5000){
                                ival = 50;
                        }else if(i < 10000){
                                ival = 100;
                        }else{
                                ival = 500;
                        }
                        ic->mapping[i] = clu_id;
                        if((i-bin_start) % ival == 0){
                                clu_id++;
                        }
                }
                /* fprintf(stderr,"%d -> %d\n", i, ic->mapping[i]); */
        }
        /* ic->clu_len = clu_id; */
        ic->clu_len = ic->mapping[len-1] + 1;


        galloc(&pd->x_axis_labels, ic->clu_len, 32);
        clu_id = 0;
        old_start = 0;
        for(int i = 0; i < len;i++){
                if(i <= bin_start){
                        snprintf(pd->x_axis_labels[clu_id],32,"%d",i);
                        old_start = i;
                        clu_id++;
                }else{

                        if(i < 2000){
                                ival = 10;
                        }else if(i < 5000){
                                ival = 50;
                        }else if(i < 10000){
                                ival = 100;
                        }else{
                                ival = 500;
                        }

                        if((i-bin_start) % ival == 0){
                                snprintf(pd->x_axis_labels[clu_id],32,"%d-%d",old_start,i);
                                old_start = i;
                                clu_id++;
                        }
                }
                /* fprintf(stderr,"%d -> %d\n", i, ic->mapping[i]); */
        }
        /* for(int i = 0; i < ic->clu_len;i++){ */
        /*         fprintf(stderr,"%s\n",pd->x_axis_labels[i]); */
        /* } */
        *interval_clu = ic;
        return OK;
ERROR:
        return FAIL;
}

int interval_clu_alloc(struct interval_clu **interval_clu, int len)
{
        struct interval_clu* ic = NULL;
        MMALLOC(ic,sizeof(struct interval_clu));

        ic->data = NULL;
        ic->n = NULL;
        ic->mapping = NULL;
        ic->L = 0;
        ic->n_clu = 0;
        ic->clu_len = 0;
        ic->org_len = len;

        galloc(&ic->mapping, ic->org_len);
        for(int i = 0; i < ic->org_len;i++){
                ic->mapping[i] = 0;
        }
        *interval_clu = ic;
        return OK;
ERROR:
        return FAIL;
}

void interval_clu_free(struct interval_clu *ic)
{
        if(ic){
                if(ic->data){
                        gfree(ic->data);
                }
                if(ic->n){
                        gfree(ic->n);
                }
                if(ic->mapping){
                        gfree(ic->mapping);
                }
                MFREE(ic);
        }
}

