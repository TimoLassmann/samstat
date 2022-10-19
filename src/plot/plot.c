#include "tld.h"

#include "plot_interval.h"

#define  PLOT_IMPORT
#include "plot.h"

struct series_data {
        char* name;
        char* label;
        char* color;
        char* type;
        int x_is_categorical;
        int vis;
        int alloc_len;
};

struct bin_data {
        double* x;
        double* y;
        double* n;
        int alloc_len;
        int len;
};



static int to_bin_data(uint64_t *x, int len,int bin_start, int bin_size, struct bin_data **bin_data);
static int bin_data_alloc(struct bin_data **bin_data, int len);
static void bin_data_free(struct bin_data *b);

#define ADD_DATA(O,S,X,Y,L) _Generic((Y),                               \
                                     uint64_t*: add_count_data,         \
                                     double*: add_real_data,            \
                                     double**: add_data_series_error    \
                )(O,S,X,Y,L)

static int write_data_normal(tld_strbuf *o, struct plot_data *d);
static int write_data_error(tld_strbuf *o, struct plot_data *d);
static int write_data_density(tld_strbuf *o, struct plot_data *d);

static int add_count_data(tld_strbuf *o,struct series_data* s_data, uint64_t *x, uint64_t* y,  int len);
static int add_real_data(tld_strbuf *o, struct series_data *s_data, double *x,
                         double *y, int len);
static int add_data_series_error(tld_strbuf *o,struct series_data* s_data, double *x, double** y,  int len);
static int fill_and_write_err(tld_strbuf *o,struct plot_data *d, int start, int stop, int plot_len);


static int add_plot_instr(tld_strbuf *o,struct plot_data *d);

static int series_data_alloc(struct series_data **series, int len);
static void series_data_free(struct series_data *s);

int plot_add(tld_strbuf *o, struct plot_data *d)
{

        char buf[256];
        for(int i= 0; i < d->L;i++){
                d->is_plot[i] = 0;
        }

        RUN(tld_append(o, "<div id=\""));
        snprintf(buf, 256,"%s_target", TLD_STR(d->id));
        RUN(tld_append(o, buf));
        RUN(tld_append(o, "\" style=\"width:100%;max-width:1400px\"></div>\n"));


        RUN(tld_append(o,"<script>\n"));
        switch(d->mod){
        case PLOT_MOD_NORMAL:
                RUN(write_data_normal(o, d));
                break;
        case PLOT_MOD_ERROR_BAR:
                RUN(write_data_error(o, d));
                break;
        case PLOT_MOD_DENSITY:
                RUN(write_data_density(o, d));
                break;
        default:
                RUN(write_data_normal(o, d));
                break;
        }

        RUN(add_plot_instr(o,d));


        RUN(tld_append(o,"</script>\n"));


        return OK;
ERROR:
        return FAIL;
}


int write_data_normal(tld_strbuf *o, struct plot_data *d)
{
        struct series_data* s_data = NULL;
        RUN(series_data_alloc(&s_data,0));
        int plot_len = 0;
        for(int i = 0; i < d->L;i++){
                for(int j = 0; j < d->len;j++){
                        if(d->data[i][j]){
                                if(j > plot_len){
                                        plot_len = j+1;
                                }
                        }
                }
        }


        for(int i = 0; i < d->L;i++){
                d->is_plot[i] = 1;
                snprintf(s_data->name,s_data->alloc_len, "%s_%d", TLD_STR(d->id),i);
                snprintf(s_data->label,s_data->alloc_len, "%s", d->series_label[i]);
                s_data->color[0] = 0;
                switch (d->type) {
                case PLOT_TYPE_SCATTER:
                        snprintf(s_data->type ,s_data->alloc_len, "scatter");
                        break;
                case PLOT_TYPE_LINES:
                        snprintf(s_data->type ,s_data->alloc_len, "lines");
                        break;
                case PLOT_TYPE_BAR:
                        snprintf(s_data->type ,s_data->alloc_len, "bar");
                        break;
                default:
                        snprintf(s_data->type ,s_data->alloc_len, "scatter");
                        break;
                }
                /* fprintf(stderr,"%s\n",s_data->type ); */
                s_data->vis = 0;
                if(d->viz == PLOT_VIZ_ALL){
                        s_data->vis = 1;
                }else{
                        if(i < d->group_size){
                                s_data->vis = 1;
                        }
                }

                /* if(i < d->group_size && d->viz == PLOT_VIZ_FIRSTGROUP){ */
                /*         s_data->vis = 1; */
                /* }else if( d->viz == PLOT_VIZ_ALL){ */
                /*         s_data->vis = 1; */
                /* }else{ */
                /*         s_data->vis = 0; */
                /* } */

                if(d->bin_size){
                        struct bin_data* b = NULL;
                        RUN(to_bin_data(d->data[i], plot_len, d->bin_start, d->bin_size, &b));
                        ADD_DATA(o, s_data, b->x, b->y ,b->len);
                        bin_data_free(b);
                }else{
                        ADD_DATA(o, s_data, NULL, d->data[i],plot_len);
                }
        }

        series_data_free(s_data);
        return OK;
ERROR:
        return FAIL;
}

int write_data_error(tld_strbuf *o, struct plot_data *d)
{
        /* struct series_data* s_data = NULL; */
        /* double* mean = NULL; */
        /* double* stderr = NULL; */
        /* double** ms = NULL; */
        /* RUN(series_data_alloc(&s_data,0)); */
        int plot_len = 0;
        for(int i = 0; i < d->L;i++){
                for(int j = 0; j < d->len;j++){
                        if(d->data[i][j]){
                                if(j > plot_len){
                                        plot_len = j+1;
                                }
                        }
                }
        }
        /* galloc(&ms, 2,plot_len); */
        /* galloc(&mean, plot_len); */
        /* galloc(&stderr, plot_len); */

        if(d->group_size < d->L){
                int start = 0;
                while(start < d->L){
                        int end = MACRO_MIN(d->L,start + d->group_size);
                        fill_and_write_err(o, d, start, end, plot_len);


                        /* RUN(tld_append(o,d->group_label[start])); */
                        start = start + d->group_size;
                }
        }else{
                fill_and_write_err(o, d, 0, d->L, plot_len);
        }
        return OK;
ERROR:
        return FAIL;
}

int write_data_density(tld_strbuf *o, struct plot_data *d)
{
        struct series_data* s_data = NULL;
        RUN(series_data_alloc(&s_data,0));
        for(int i = 0; i < d->L;i++){
                double** density = NULL;
                double* x = NULL;
                double* y = NULL;
                galloc(&x,MACRO_MAX(10,d->len));
                galloc(&y,MACRO_MAX(10,d->len));
                int c = 0;
                for(int j = 0; j < d->len;j++){
                        if(d->data[i][j]){
                                /* LOG_MSG("Adding %d",i); */
                                x[c] = (double) j;
                                y[c] = (double) d->data[i][j];
                                /* LOG_MSG("%f %f", x[c],y[c]); */
                                c++;
                        }
                        /* dat[i] = (double) l->data[i]; */
                }

                if(c == 0){
                        ERROR_MSG("Not a single sequence found in length distribution");
                }else if(c == 1){
                        galloc(&density, 1,2);
                        density[0][0] = x[0];
                        density[0][1] = 1.0;
                }else {
                        tld_kde_count_pdf(x,y,c,  &density);
                }

                for(int j = 0; j < c;j++){
                        x[j] = density[j][0];
                        y[j] = density[j][1];
                }
                snprintf(s_data->name,s_data->alloc_len, "%s_%d", TLD_STR(d->id),i);
                snprintf(s_data->label,s_data->alloc_len, "%s", d->series_label[i]);
                s_data->color[0] = 0;
                snprintf(s_data->type ,s_data->alloc_len, "lines");
                /* fprintf(stderr,"%s\n",s_data->type ); */
                s_data->vis = 1;

                if(d->bin_size){
                        struct bin_data* b = NULL;
                        RUN(to_bin_data(d->data[i], c, d->bin_start, d->bin_size, &b));
                        ADD_DATA(o, s_data, b->x, b->y ,b->len);
                        bin_data_free(b);
                }else{
                        ADD_DATA(o, s_data, x, y,c);
                }

                d->is_plot[i] = 1;
                gfree(x);
                gfree(y);
                gfree(density);
        }

        series_data_free(s_data);
        return OK;
ERROR:
        series_data_free(s_data);
        return FAIL;

}

int fill_and_write_err(tld_strbuf *o,struct plot_data *d, int start, int stop, int plot_len)
{
        struct series_data* s_data = NULL;
        double **ms = NULL;
        RUN(series_data_alloc(&s_data,0));
        galloc(&ms, plot_len,2);
        for(int i = 0; i < plot_len;i++){
                ms[i][0] = 0.0;
                ms[i][1] = 0.0;
        }
        for(int i = 0 ; i < plot_len;i++){
                double total = 0.0;
                double n = 0.0;
                for(int j= start; j < stop;j++){
                        total += (j-start) * d->data[j][i];
                        n+= d->data[j][i];
                }
                ms[i][0] = total / n;
                total = 0.0;
                n = 0.0;
                for(int j = start; j < d->L;j++){
                        double s =  ms[i][0] - (double) (j-start) ;
                        total += d->data[j][i] *  s * s;
                        n += d->data[j][i];
                }
                ms[i][1] = sqrt(total / n) / sqrt(n);
        }
        d->is_plot[start] = 1;
        snprintf(s_data->name,s_data->alloc_len, "%s_%d", TLD_STR(d->id),start);
        snprintf(s_data->label,s_data->alloc_len, "%s", d->group_label[start]);
        s_data->color[0] = 0;
        switch (d->type) {
        case PLOT_TYPE_SCATTER:
                snprintf(s_data->type ,s_data->alloc_len, "scatter");
                break;
        case PLOT_TYPE_LINES:
                snprintf(s_data->type ,s_data->alloc_len, "lines");
                break;
        case PLOT_TYPE_BAR:
                snprintf(s_data->type ,s_data->alloc_len, "bar");
                break;
        default:
                snprintf(s_data->type ,s_data->alloc_len, "scatter");
                break;
        }

        s_data->vis = 1;
        /* if(start == 0){ */
        /*         s_data->vis = 1; */
        /* } */
        ADD_DATA(o, s_data, NULL, ms,plot_len);
        series_data_free(s_data);
        gfree(ms);
        return OK;
ERROR:
        series_data_free(s_data);
        gfree(ms);
        return FAIL;
}

int add_plot_instr(tld_strbuf *o, struct plot_data *d)
{
        char buf[256];
        /* int viz = d->type >> PLOT_TYPE_SHIFT; */

        /* Data */
        RUN(tld_append(o,"var "));
        snprintf(buf, 256,"%s_data", TLD_STR(d->id));
        RUN(tld_append(o,buf ));
        RUN(tld_append(o," = ["));
        for(int i = 0; i < d->L;i++){
                if(d->is_plot[i]){
                        snprintf(buf, 256,"%s_%d,", TLD_STR(d->id),i);
                        RUN(tld_append(o,buf));
                }

                /* snprintf(s_data->name,s_data->alloc_len, "%s_%d,", TLD_STR(d->id),i); */

        }
        o->len--;
        RUN(tld_append(o,"];\n"));
        /* Layout */
        RUN(tld_append(o,"var "));
        snprintf(buf, 256,"%s_layout", TLD_STR(d->id));
        RUN(tld_append(o,buf ));
        RUN(tld_append(o," = {\n"));
        RUN(tld_append(o,"boxmode: 'group',\n"));
        if(d->mod != PLOT_MOD_ERROR_BAR){
                RUN(tld_append(o,"barmode: 'stack',\n"));
        }
        RUN(tld_append(o,"title: '"));
        RUN(tld_append(o,TLD_STR(d->title)));
        RUN(tld_append(o,"',\n"));
        /* Buttom?  */

        if(d->group_size < d->L && d->viz == PLOT_VIZ_FIRSTGROUP){
                RUN(tld_append(o,"updatemenus: [{\n"));
                RUN(tld_append(o,"x: 0.5,\n"));
                RUN(tld_append(o,"xanchor: 'center',\n"));
                RUN(tld_append(o,"y: 1.09,\n"));
                RUN(tld_append(o,"yanchor: 'top',\n"));
                RUN(tld_append(o,"buttons: [\n"));
                int start = 0;
                while(start < d->L){
                        if(start){
                                RUN(tld_append(o,",\n"));
                        }
                        RUN(tld_append(o,"{\n"));
                        RUN(tld_append(o,"method: 'restyle',\n"));
                        RUN(tld_append(o,"args: ['visible', [\n"));
                        for(int i = 0; i < d->L;i++){
                                if(i >= start && i < start+d->group_size){
                                        RUN(tld_append(o,"true,"));
                                }else{
                                        RUN(tld_append(o,"false,"));
                                }
                        }
                        o->len--;
                        RUN(tld_append(o,"]],\n"));
                        RUN(tld_append(o,"label: '"));
                        RUN(tld_append(o,d->group_label[start]));
                        RUN(tld_append(o,"'\n"));

                        RUN(tld_append(o,"}\n"));
                        start = start + d->group_size;
                }
                RUN(tld_append(o,",]}],\n"));
        }


        RUN(tld_append(o,"xaxis: {\n"));
        /* RUN(tld_append(o,"type: 'category',\n")); */
        if(d->len > 100){
                RUN(tld_append(o,"range: [0, 75],\n"));
        }
        snprintf(buf, 256,"title: '%s'\n", TLD_STR(d->xlabel));

        RUN(tld_append(o,buf ));

        if(d->len > 100){
                RUN(tld_append(o,",\n"));

                RUN(tld_append(o,"rangeslider: {}\n"));
        }
        RUN(tld_append(o,"},\n"));
        RUN(tld_append(o,"yaxis: {\n"));
        RUN(tld_append(o,"rangemode: \"nonnegative\",\n"));
        snprintf(buf, 256,"title: '%s'\n", TLD_STR(d->ylabel));
        RUN(tld_append(o,buf ));
        RUN(tld_append(o,"}\n"));
        RUN(tld_append(o,"};\n"));

        RUN(tld_append(o,"Plotly.newPlot("));
        snprintf(buf, 256,"%s_target", TLD_STR(d->id));
        RUN(tld_append(o,buf)) ;
        RUN(tld_append(o,","));
        snprintf(buf, 256,"%s_data", TLD_STR(d->id));
        RUN(tld_append(o,buf)) ;
        RUN(tld_append(o,","));
        snprintf(buf, 256,"%s_layout", TLD_STR(d->id));
        RUN(tld_append(o,buf ));
        RUN(tld_append(o,",\n"));
        RUN(tld_append(o,"{toImageButtonOptions : {\n"));
        RUN(tld_append(o,"format : 'svg',\n"));
        RUN(tld_append(o,"filename : '"));
        RUN(tld_append(o,TLD_STR(d->title)));
        RUN(tld_append(o,"',\n"));


        RUN(tld_append(o,"height : 500,\n"));
        RUN(tld_append(o,"width : 700,\n"));
        RUN(tld_append(o,"scale : 1}}\n"));


        RUN(tld_append(o,")\n"));




/*         var base_comp_layout = { */
/* 138-grid: {rows: 1, columns: 2, pattern: 'independent'}, */
/* 139-barmode: 'stack'}; */
/* 140-Plotly.newPlot('basecomp', base_comp_data,{barmode: 'stack', */
/* 141-title: 'Base Composition', */
/* 142-xaxis: { */
/* 143-title: 'Length' */
/* 144-}, */
/* 145-yaxis: { */
/* 146-title: 'Counts' */
/* 147-}, */
/* -- */
/* 182:var quallayout = {boxmode: 'group', */
/* 183-xaxis: { */
/* 184-title: 'Length' */
/* 185-}, */
/* 186-yaxis: { */
/* 187-title: 'Base Quality' */
/* 188-} */
/* 189-}; */


/*         RUN(tld_append(o,"Plotly.newPlot('")); */

              /* tld_strbuf_free_space(tld_strbuf *b) */
/*         var qualdata = [qualtrace0,qualtrace2,qualtrace3]; */
/* var quallayout = {boxmode: 'group', */
/* xaxis: { */
/* title: 'Length' */
/* }, */
/* yaxis: { */
/* title: 'Base Quality' */
/* } */
/* }; */
/* Plotly.newPlot('qualcomp', qualdata, quallayout, { */
/*   toImageButtonOptions : { */
/*     format : 'svg', */
/*     filename : 'ENCFF039WAK.bam_BaseQuality', */
/*     height : 500, */
/*     width : 700, */
/*     scale : 1 */
/*   } */
/* }); */
        return OK;
ERROR:
        return FAIL;
}

int plot_data_alloc(struct plot_data** plot_data,int x, int y)
{

        struct plot_data* pd = NULL;
        MMALLOC(pd, sizeof(struct plot_data));
        pd->data = NULL;
        pd->series_label = NULL;
        pd->group_label = NULL;
        pd->is_plot = NULL;
        pd->len = x;
        pd->L = y;
        pd->title = NULL;
        pd->id = NULL;
        pd->xlabel = NULL;
        pd->ylabel = NULL;
        pd->save_file_name = NULL;
        /* pd->type = 0; */
        pd->viz = pd->L;
        pd->type = PLOT_TYPE_SCATTER;
        pd->mod = PLOT_MOD_NORMAL;
        pd->bin_size = 0;
        pd->bin_start= 0;

        RUN(tld_strbuf_alloc(&pd->title, 256));
        RUN(tld_strbuf_alloc(&pd->id, 256));
        RUN(tld_strbuf_alloc(&pd->xlabel, 256));
        RUN(tld_strbuf_alloc(&pd->ylabel, 256));
        RUN(tld_strbuf_alloc(&pd->save_file_name, 256));

        galloc(&pd->data,pd->L,pd->len );
        galloc(&pd->is_plot,pd->L);
        galloc(&pd->series_label,pd->L,256);
        galloc(&pd->group_label,pd->L,256);

        for(int i =0;i < pd->L;i++){
                pd->is_plot[i] = 0;
                for(int j = 0; j < pd->len;j++){
                        pd->data[i][j] = 0;
                }
        }

        *plot_data = pd;
        return OK;
ERROR:
        plot_data_free(pd);
        return FAIL;
}

int plot_data_resize_len(struct plot_data* pd,int x)
{
        if(pd->len > x){
                ERROR_MSG("New len is smaller than before!");
        }
        uint64_t** new = NULL;
        int old_len = pd->len;
        pd->len = x;
        galloc(&new,pd->L,pd->len);

        for(int i = 0; i < pd->L;i++){
                for(int j = 0; j < old_len;j++){
                        new[i][j] = pd->data[i][j];
                }
                for(int j = old_len; j < pd->len;j++){
                        new[i][j] = 0;
                }

        }

        gfree(pd->data);
        pd->data = new;
        return OK;
ERROR:
        return FAIL;
}


void plot_data_free(struct plot_data *pd)
{
        if(pd){
                tld_strbuf_free(pd->title);
                tld_strbuf_free(pd->id);
                tld_strbuf_free(pd->xlabel);
                tld_strbuf_free(pd->ylabel);
                tld_strbuf_free(pd->save_file_name);

                gfree(pd->data);
                gfree(pd->series_label);
                gfree(pd->is_plot);
                gfree(pd->group_label);
                MFREE(pd);
        }

}

int series_data_alloc(struct series_data **series, int len)
{
        struct series_data* s = NULL;
        MMALLOC(s, sizeof(struct series_data));
        s->name = NULL;
        s->label = NULL;
        s->color = NULL;
        s->type = NULL;
        s->vis = 0;
        s->x_is_categorical = 0;
        if(len){
                s->alloc_len = len;
        }else{
                s->alloc_len = 256;
        }
        galloc(&s->name,s->alloc_len);
        galloc(&s->label,s->alloc_len);
        galloc(&s->color,s->alloc_len);
        galloc(&s->type,s->alloc_len);
        *series = s;
        return OK;
ERROR:
        series_data_free(s);
        return FAIL;
}

void series_data_free(struct series_data *s)
{
        if(s){
                if(s->name){
                        gfree(s->name);
                }
                if(s->label){
                        gfree(s->label);
                }
                if(s->color){
                        gfree(s->color);
                }
                if(s->type){
                        gfree(s->type);
                }
                MFREE(s);
        }
}


int add_count_data(tld_strbuf *o,struct series_data* s_data , uint64_t *x, uint64_t* y,  int len)
{
        ASSERT(y != NULL, "no y");
        char buf[256];
        snprintf(buf, 256,"var %s  = {\n",s_data->name);
        RUN(tld_append(o,buf));
        if(x == NULL){
                RUN(tld_append(o,"x: ["));
                snprintf(buf, 256,"%d",1);
                RUN(tld_append(o,buf));
                for(int i = 1 ; i < len;i++){
                        snprintf(buf, 256,",%d",i+1);
                        RUN(tld_append(o,buf));
                }
                RUN(tld_append(o,"],\n"));
        }else{
                RUN(tld_append(o,"x: ["));

                for(int i = 0 ; i < len;i++){
                        snprintf(buf, 256,"%"PRId64",",x[i]);
                        RUN(tld_append(o,buf));
                }
                o->len--;
                RUN(tld_append(o,"],\n"));
        }
        RUN(tld_append(o,"y: ["));
        snprintf(buf, 256,"%"PRId64"", y[0]);
        RUN(tld_append(o,buf));
        for(int i = 1 ; i < len;i++){
                snprintf(buf, 256,",%"PRId64"",y[i]);
                RUN(tld_append(o,buf));
        }
        RUN(tld_append(o,"],\n"));
        if(s_data->label){
                snprintf(buf, 256,"name: '%s',\n", s_data->label);
                RUN(tld_append(o,buf));
        }
        if(s_data->color[0] != 0){
                /* snprintf(buf, 256,"name: %s,\n", label); */
                /* RUN(tld_append(o,buf)); */
        }
        if(s_data->type){
                snprintf(buf, 256,"type: '%s',\n", s_data->type);
                RUN(tld_append(o,buf));
        }else{
                RUN(tld_append(o,"type: 'scatter',\n"));
        }
        if(!s_data->vis){
                RUN(tld_append(o,"visible: false\n"));
        }else{
                RUN(tld_append(o,"visible: true\n"));
        }
        RUN(tld_append(o,"};\n"));
        return OK;
ERROR:
        return FAIL;
}

int add_real_data(tld_strbuf *o,struct series_data* s_data, double *x, double* y,  int len)
{
        ASSERT(y != NULL, "no y");
        char buf[256];
        snprintf(buf, 256,"var %s  = {\n",s_data->name);
        RUN(tld_append(o,buf));
        if(x == NULL){
                RUN(tld_append(o,"x: ["));
                snprintf(buf, 256,"%d",1);
                RUN(tld_append(o,buf));
                for(int i = 1 ; i < len;i++){
                        snprintf(buf, 256,",%d",i+1);
                        RUN(tld_append(o,buf));
                }
                RUN(tld_append(o,"],\n"));
        }else{
                RUN(tld_append(o,"x: ["));

                for(int i = 0 ; i < len;i++){
                        snprintf(buf, 256,"\"%f\",",x[i]);
                        RUN(tld_append(o,buf));
                }
                o->len--;
                RUN(tld_append(o,"],\n"));
        }
        RUN(tld_append(o,"y: ["));
        /* snprintf(buf, 256,"%f", y[0]); */
        /* RUN(tld_append(o,buf)); */
        for(int i = 0 ; i < len;i++){
                snprintf(buf, 256,"%f,",y[i]);
                RUN(tld_append(o,buf));
        }
        o->len--;
        RUN(tld_append(o,"],\n"));
        if(s_data->label){
                snprintf(buf, 256,"name: '%s',\n", s_data->label);
                RUN(tld_append(o,buf));
        }
        if(s_data->color[0] != 0){
                /* snprintf(buf, 256,"name: %s,\n", label); */
                /* RUN(tld_append(o,buf)); */
        }
        if(s_data->type){
                snprintf(buf, 256,"type: '%s',\n", s_data->type);
                RUN(tld_append(o,buf));
        }else{
                RUN(tld_append(o,"type: 'scatter',\n"));
        }
        if(!s_data->vis){
                RUN(tld_append(o,"visible: false\n"));
        }else{
                RUN(tld_append(o,"visible: true\n"));
        }

        RUN(tld_append(o,"};\n"));
        return OK;
ERROR:
        return FAIL;
}


int add_data_series_error(tld_strbuf *o,struct series_data* s_data, double *x, double** y,  int len)
{
        ASSERT(y != NULL, "no y");
        char buf[256];
        snprintf(buf, 256,"var %s  = {\n",s_data->name);
        RUN(tld_append(o,buf));
        if(x == NULL){
                RUN(tld_append(o,"x: ["));
                snprintf(buf, 256,"%d",1);
                RUN(tld_append(o,buf));
                for(int i = 1 ; i < len;i++){
                        snprintf(buf, 256,",%d",i+1);
                        RUN(tld_append(o,buf));
                }
                RUN(tld_append(o,"],\n"));
        }else{
                RUN(tld_append(o,"x: ["));

                for(int i = 0 ; i < len;i++){
                        snprintf(buf, 256,"%f,",x[i]);
                        RUN(tld_append(o,buf));
                }
                o->len--;
                RUN(tld_append(o,"],\n"));
        }
        RUN(tld_append(o,"y: ["));
        /* snprintf(buf, 256,"%f", y[0]); */
        /* RUN(tld_append(o,buf)); */
        for(int i = 0 ; i < len;i++){
                snprintf(buf, 256,"%f,",y[i][0]);
                RUN(tld_append(o,buf));
        }
        o->len--;
        RUN(tld_append(o,"],\n"));

        RUN(tld_append(o,"error_y: {\n"));
        RUN(tld_append(o,"type: 'data',\n"));
        RUN(tld_append(o,"array : ["));
        for(int i = 0 ; i < len;i++){
                snprintf(buf, 256,"%f,",y[i][1]);
                RUN(tld_append(o,buf));
        }
        o->len--;
        RUN(tld_append(o,"],\n"));
        RUN(tld_append(o,"visible: true\n"));
        RUN(tld_append(o,"},\n"));


        if(s_data->label){
                snprintf(buf, 256,"name: '%s',\n", s_data->label);
                RUN(tld_append(o,buf));
        }
        if(s_data->color[0] != 0){
                /* snprintf(buf, 256,"name: %s,\n", label); */
                /* RUN(tld_append(o,buf)); */
        }
        if(s_data->type){
                snprintf(buf, 256,"type: '%s',\n", s_data->type);
                RUN(tld_append(o,buf));
        }else{
                RUN(tld_append(o,"type: 'scatter',\n"));
        }
        if(!s_data->vis){
                RUN(tld_append(o,"visible: false\n"));
        }else{
                RUN(tld_append(o,"visible: true\n"));
        }

        RUN(tld_append(o,"};\n"));
        return OK;
ERROR:
        return FAIL;
}

int html_header(tld_strbuf *out_buffer, char *filename)
{
        RUN(tld_append(out_buffer, "<!DOCTYPE html>\n"));
        RUN(tld_append(out_buffer, "<html>\n"));

        RUN(tld_append(out_buffer, "<head>\n"));

        RUN(tld_append(out_buffer, "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"/>\n"))
        RUN(tld_append(out_buffer, "<script src=\"https://d3js.org/d3.v4.js\"></script>\n"));
        RUN(tld_append(out_buffer, "<script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>\n"));


        if(filename){
                RUN(tld_append(out_buffer, "<title>\n"));
                RUN(tld_append(out_buffer, filename));
                RUN(tld_append(out_buffer, "</title>\n"));
        }

        RUN(tld_append(out_buffer, "</head>\n"));

        /* RUN(tld_append(out_buffer, "<title>%s</title>\n", "ARFFF"); */

        RUN(tld_append(out_buffer, "<style>\n"));



        //RUN(tld_append(out_buffer, "canvas{\n"));
        //RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "body {\n"));
        RUN(tld_append(out_buffer, "	margin: 0 auto;\n"));
        RUN(tld_append(out_buffer, "	padding: 22px 0;\n"));
        RUN(tld_append(out_buffer, "	font-family: 'Source Sans Pro', sans-serif;\n"));
        RUN(tld_append(out_buffer, "	font-weight: normal;\n"));
        RUN(tld_append(out_buffer, "	width: 940px;\n"));
        //RUN(tld_append(out_buffer, "font-family: \"Open Sans\";\n"));
        //RUN(tld_append(out_buffer, "font: 14px Arial,Helvetica, sans-serif;\n"));
        RUN(tld_append(out_buffer, "	background: #F0F0F0;\n"));
        RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "h2 {\n"));
        RUN(tld_append(out_buffer, "	font-weight: 200;\n"));
        RUN(tld_append(out_buffer, "	font-size: 1.5em;\n"));
        RUN(tld_append(out_buffer, "	border-bottom: solid 1px black;\n"));
        RUN(tld_append(out_buffer, "	padding-top: 35px;\n"));
        RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "p {\n"));
        RUN(tld_append(out_buffer, "	line-height:1.4;\n"));
        RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "footer {\n"));
        RUN(tld_append(out_buffer, "	position: absolute;\n"));
        RUN(tld_append(out_buffer, "	left: 0;\n"));
        RUN(tld_append(out_buffer, "	width: 100%;\n"));
        RUN(tld_append(out_buffer, "	background: #222;\n"));
        RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "footer div {\n"));
        RUN(tld_append(out_buffer, "	display: table;\n"));
        RUN(tld_append(out_buffer, "	margin: 0 auto;\n"));
        RUN(tld_append(out_buffer, "	padding: 44px 0;\n"));
        RUN(tld_append(out_buffer, "	width: 940px;\n"));
        RUN(tld_append(out_buffer, "	color: #777;\n"));
        RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "footer div section {\n"));
        RUN(tld_append(out_buffer, "	display: table-cell;\n"));
        RUN(tld_append(out_buffer, "	width: 300px;\n"));
        RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "footer div #about, footer div #blogroll {\n"));
        RUN(tld_append(out_buffer, "	padding-right: 20px;\n"));
        RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "footer h3 {\n"));
        RUN(tld_append(out_buffer, "	color: #FFF;\n"));
        RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "footer a {\n"));
        RUN(tld_append(out_buffer, "	color: #999;\n"));
        RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "footer a:hover {\n"));
        RUN(tld_append(out_buffer, "	color: #FFF;\n"));
        RUN(tld_append(out_buffer, "	text-decoration: none;\n"));
        RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "footer ul {\n"));
        RUN(tld_append(out_buffer, "    margin: 0 0 0 40px;\n"));
        RUN(tld_append(out_buffer, "	list-style: square;\n"));
        RUN(tld_append(out_buffer, "    color: #565656;\n"));
        RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "footer ul li a {\n"));
        RUN(tld_append(out_buffer, "    display: block;\n"));
        RUN(tld_append(out_buffer, "}\n"));

        /* RUN(tld_append(out_buffer, "canvas  {\n")); */
        /* RUN(tld_append(out_buffer, "	padding-left: 0;\n")); */
        /* RUN(tld_append(out_buffer, "	padding-right: 0;\n")); */
        /* RUN(tld_append(out_buffer, "	margin-left: 0;\n")); */
        /* RUN(tld_append(out_buffer, "	margin-right: auto;\n")); */
        /* RUN(tld_append(out_buffer, "	display: block;\n")); */
        /* RUN(tld_append(out_buffer, "	float: left;\n")); */
        /* RUN(tld_append(out_buffer, "}\n")); */


        RUN(tld_append(out_buffer, "nav {\n"));
        //RUN(tld_append(out_buffer, "position: absolute;\n"));
        RUN(tld_append(out_buffer, "left: 0;\n"));
        RUN(tld_append(out_buffer, "width: 100%%;\n"));
        RUN(tld_append(out_buffer, "background: #75787B;\n"));
        RUN(tld_append(out_buffer, "color: #FFF;\n"));
        RUN(tld_append(out_buffer, "}\n"));




        RUN(tld_append(out_buffer, "</style>\n"));
        RUN(tld_append(out_buffer, "<body>\n"));

        return OK;
ERROR:
        return FAIL;
}

int html_end(tld_strbuf *out_buffer)
{
        RUN(tld_append(out_buffer, "</body>\n"));
        RUN(tld_append(out_buffer, "</html>\n"));
        return OK;
ERROR:
        return FAIL;
}

int to_bin_data(uint64_t *x, int len,int bin_start, int bin_size, struct bin_data **bin_data)
{
        struct bin_data* b = NULL;
        int idx;
        int max_idx = 0;
        int interval = 0;
        get_plot_interval(len, bin_start, 50, &interval);

        fprintf(stderr,"interval: %d\n", interval);

        RUN(bin_data_alloc(&b, len));
        idx = 0;
        for(int i = 0; i < len;i++){

                if(i <= bin_start){
                        b->x[idx] = i;
                        b->y[idx] = x[i];
                        b->n[idx] = 1.0;
                        fprintf(stderr,"%d -> %d single \n", i, idx);
                        idx++;
                }else{
                        idx = (i-1) / bin_size + bin_start +1;
                        fprintf(stderr,"%d -> %d\n", i, idx);
                        b->x[idx] = MACRO_MAX(i, b->x[idx]);
                        b->y[idx] += x[i];
                        b->n[idx] += 1.0;

                }
                max_idx = MACRO_MAX(max_idx, idx);
        }
        max_idx+= 1;
        fprintf(stderr,"Max idx : %d\n", max_idx);
        b->len = max_idx;
        for(int i = 0; i < b->len;i++){
                ASSERT(b->n[i] > 0.0,"N can never be zero!");
                fprintf(stderr,"%d:x:%f %f %f -> %f \n", i,b->x[i], b->y[i] , b->n[i], b->y[i] / b->n[i]);
                b->y[i] = b->y[i] / b->n[i];
        }
        exit(0);
        *bin_data = b;
        return OK;
ERROR:
        bin_data_free(b);
        return FAIL;
}

int bin_data_alloc(struct bin_data **bin_data, int len)
{
        struct bin_data* b = NULL;
        MMALLOC(b, sizeof(struct bin_data));
        b->len = 0;
        b->alloc_len = len;
        b->x = NULL;
        b->y = NULL;
        b->n = NULL;
        galloc(&b->x, b->alloc_len);
        galloc(&b->y, b->alloc_len);
        galloc(&b->n, b->alloc_len);
        for(int i = 0; i < b->alloc_len;i++){
                b->x[i] = 0;
                b->y[i] = 0;
                b->n[i] = 0;
        }
        *bin_data = b;

        return OK;
ERROR:
        bin_data_free(b);
        return FAIL;
}

void bin_data_free(struct bin_data *b)
{
        if(b){
                if(b->x){
                        gfree(b->x);
                }
                if(b->y){
                        gfree(b->y);
                }
                if(b->n){
                        gfree(b->n);
                }
                MFREE(b);
        }

}
