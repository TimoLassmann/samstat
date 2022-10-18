#include "tld.h"



#define  PLOT_IMPORT
#include "plot.h"

struct series_data {
        char* name;
        char* label;
        char* color;
        char* type;
        int vis;
        int alloc_len;
};

#define ADD_DATA(O,S,X,Y,L) _Generic((Y),                       \
                                     uint64_t*: add_count_data, \
                                     double*: add_real_data     \
                )(O,S,X,Y,L)

static int add_count_data(tld_strbuf *o,struct series_data* s_data, uint64_t *x, uint64_t* y,  int len);
static int add_real_data(tld_strbuf *o,struct series_data* s_data, double *x, double* y,  int len);

static int add_plot_instr(tld_strbuf *o,struct plot_data *d);

static int series_data_alloc(struct series_data **series, int len);
static void series_data_free(struct series_data *s);

int plot_add(tld_strbuf *o, struct plot_data *d)
{
        struct series_data* s_data = NULL;
        char buf[256];
        RUN(series_data_alloc(&s_data,0));


        /* RUN(tld_append(o, "<div id=\"")); */
        /* snprintf(buf, 256,"delerrcomp"); */
        /* RUN(tld_append(o, TLD_STR(d->id))); */
        /* RUN(tld_append(o, "\" style=\"width:100%;max-width:1400px\"></div>\n")); */
        /* RUN(tld_append(o,"<script>\n")); */

        int plot_type = d->type & PLOT_TYPE_MASK;
        int viz = d->type >> PLOT_TYPE_SHIFT;

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


        RUN(tld_append(o, "<div id=\""));
        snprintf(buf, 256,"%s_target", TLD_STR(d->id));
        RUN(tld_append(o, buf));
        RUN(tld_append(o, "\" style=\"width:100%;max-width:1400px\"></div>\n"));


        RUN(tld_append(o,"<script>\n"));
        for(int i = 0; i < d->L;i++){
                snprintf(s_data->name,s_data->alloc_len, "%s_%d", TLD_STR(d->id),i);
                snprintf(s_data->label,s_data->alloc_len, "%s", d->series_label[i]);
                s_data->color[0] = 0;
                switch (plot_type) {
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
                fprintf(stderr,"%s\n",s_data->type );
                s_data->vis = 0;
                if(i < viz){
                        s_data->vis = 1;
                }
                ADD_DATA(o, s_data, NULL, d->data[i],plot_len);
        }

        RUN(add_plot_instr(o,d));


        RUN(tld_append(o,"</script>\n"));



        series_data_free(s_data);

        return OK;
ERROR:
        return FAIL;
}

int add_plot_instr(tld_strbuf *o, struct plot_data *d)
{
        char buf[256];
        int viz = d->type >> PLOT_TYPE_SHIFT;
        /* Data */
        RUN(tld_append(o,"var "));
        snprintf(buf, 256,"%s_data", TLD_STR(d->id));
        RUN(tld_append(o,buf ));
        RUN(tld_append(o," = ["));
        for(int i = 0; i < d->L;i++){
                snprintf(buf, 256,"%s_%d,", TLD_STR(d->id),i);
                RUN(tld_append(o,buf));

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
        RUN(tld_append(o,"barmode: 'stack',\n"));
        RUN(tld_append(o,"title: '"));
        RUN(tld_append(o,TLD_STR(d->title)));
        RUN(tld_append(o,"',\n"));
        /* Buttom?  */

        if(viz < d->L){
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
                                if(i >= start && i < start+viz){
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
                        start = start + viz;
                }
                RUN(tld_append(o,",]}],\n"));
        }


        RUN(tld_append(o,"xaxis: {\n"));
        snprintf(buf, 256,"title: '%s'\n", TLD_STR(d->xlabel));
        RUN(tld_append(o,buf ));
        RUN(tld_append(o,"},\n"));
        RUN(tld_append(o,"yaxis: {\n"));
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
        pd->len = x;
        pd->L = y;
        pd->title = NULL;
        pd->id = NULL;
        pd->xlabel = NULL;
        pd->ylabel = NULL;
        pd->save_file_name = NULL;
        pd->type = 0;

        RUN(tld_strbuf_alloc(&pd->title, 256));
        RUN(tld_strbuf_alloc(&pd->id, 256));
        RUN(tld_strbuf_alloc(&pd->xlabel, 256));
        RUN(tld_strbuf_alloc(&pd->ylabel, 256));
        RUN(tld_strbuf_alloc(&pd->save_file_name, 256));

        galloc(&pd->data,pd->L,pd->len );
        galloc(&pd->series_label,pd->L,256);
        galloc(&pd->group_label,pd->L,256);
        *plot_data = pd;
        return OK;
ERROR:
        plot_data_free(pd);
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
