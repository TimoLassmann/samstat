#include "core/tld-core.h"
#include "seq/tld-seq.h"
#include "stats/basic.h"
#include "string/str.h"
#include "tld.h"
#include <math.h>
#include <stdio.h>

#define REPORT_IMPORT
#include "report.h"

static int report_header(tld_strbuf *out_buffer);
static int report_footer(tld_strbuf *out_buffer);

static int mapping_quality_overview_section(tld_strbuf *o, struct metrics *m);
static int base_composition_section(tld_strbuf *o, struct metrics *m);
static int base_quality_section(tld_strbuf *o, struct metrics *m);

int report_header(tld_strbuf *out_buffer)
{
        RUN(tld_append(out_buffer, "<!DOCTYPE html>\n"));
        RUN(tld_append(out_buffer, "<html>\n"));
        RUN(tld_append(out_buffer, "<script src=\"https://d3js.org/d3.v4.js\"></script>\n"));
        RUN(tld_append(out_buffer, "<script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>\n"));
        RUN(tld_append(out_buffer, "<body>\n"));
        return OK;
ERROR:
        return FAIL;
}

int report_footer(tld_strbuf *out_buffer)
{
        RUN(tld_append(out_buffer, "</body>\n"));
        RUN(tld_append(out_buffer, "</html>\n"));
        return OK;
ERROR:
        return FAIL;
}

int mapping_quality_overview_section(tld_strbuf *o, struct metrics *m)
{
        char buf[256];
        if(!m->is_aligned){
                return OK;
        }

        RUN(tld_append(o, "<h2>Mapping Stats</h2>\n"));
        RUN(tld_append(o, "<div id=\"MappingStats\" style=\"width:100%;max-width:1400px\"></div>\n"));

        RUN(tld_append(o,"<script>\n"));
        RUN(tld_append(o, "var mapping_data = [{ \n"));
        RUN(tld_append(o, "type: \"pie\",\n"));
        RUN(tld_append(o, "values: ["));
        snprintf(buf, 256, "%d",m->seq_comp_R1[0]->n_counts);
        RUN(tld_append(o,buf));
        for(int i = 1; i < 6;i++){
                snprintf(buf, 256, ",%d",m->seq_comp_R1[i]->n_counts);
                RUN(tld_append(o,buf));
        }
        RUN(tld_append(o, "],"));
        RUN(tld_append(o, "labels: ["));
        snprintf(buf, 256, "\"%s (%d sequences)\"",m->mapq_map->description[0],m->seq_comp_R1[0]->n_counts);
        RUN(tld_append(o,buf));
        for(int i = 1; i < m->mapq_map->n_bin;i++){
                snprintf(buf, 256, ",\"%s (%d sequences)\"",m->mapq_map->description[i],m->seq_comp_R1[i]->n_counts );
                RUN(tld_append(o,buf));
        }
        RUN(tld_append(o, "],\n"));

        RUN(tld_append(o, "textinfo: \"label+percent\",\n"));
        RUN(tld_append(o, "textposition: \"outside\",\n"));
        RUN(tld_append(o, "automargin: true\n"));

        RUN(tld_append(o, "}]\n"));

        RUN(tld_append(o, "var mapping_layout = {\n"));
        RUN(tld_append(o, "height: 400,\n"));
        RUN(tld_append(o, "width: 400,\n"));
        RUN(tld_append(o, "margin: {\"t\": 0, \"b\": 0, \"l\": 0, \"r\": 0},\n"));
        RUN(tld_append(o, "showlegend: true\n"));
        RUN(tld_append(o, "}\n"));

        RUN(tld_append(o, "Plotly.newPlot('MappingStats', mapping_data, mapping_layout)\n"));

        RUN(tld_append(o,"</script>\n"));
        return OK;
ERROR:
        return FAIL;
}

int base_quality_section(tld_strbuf *o, struct metrics *m)
{
        struct qual_composition* q = NULL;
        char buf[256];
        double* mean = NULL;
        double* stderr = NULL;

        q = m->qual_comp_R1[0];

        galloc(&mean, q->len);
        galloc(&stderr, q->len);

        RUN(tld_append(o, "<h2>Base quality distribution</h2>\n"));
        RUN(tld_append(o, "<div id=\"qualcomp\" style=\"width:100%;max-width:1400px\"></div>\n"));
        RUN(tld_append(o, "<script>\n"));
        RUN(tld_append(o, "var basex = [\n"));
        for(int i = 0 ; i < q->len;i++){
                for(int j = 0; j < q->L;j++){
                        snprintf(buf, 256,"%d,",i);
                        RUN(tld_append(o,buf));
                }
        }
        o->len--;
        RUN(tld_append(o, "]\n"));

        for(int mapq_idx = 0; mapq_idx < 6; mapq_idx++){
                q = m->qual_comp_R1[mapq_idx];
                if(q->n_counts > 0){

                        /* calc mean */

                        for(int i = 0 ; i < q->len;i++){
                                double total = 0.0;

                                double n = 0.0;
                                for(int j = 0; j < q->L;j++){
                                        total += j * q->data[i][j];
                                        n+= q->data[i][j];
                                }
                                mean[i] = total / n;
                                total = 0.0;
                                n = 0.0;
                                for(int j = 0; j < q->L;j++){
                                        double s = mean[i] - (double) j ;
                                        total += q->data[i][j] *  s * s;
                                        n+= q->data[i][j];
                                }
                                stderr[i] = sqrt(total / n) / sqrt(n);
                                LOG_MSG("%d %f %f", i, mean[i], stderr[i]);
                        }

                        snprintf(buf, 256,"var qualtrace%d = {\n",mapq_idx);
                        RUN(tld_append(o,buf));

                        RUN(tld_append(o,"name: '"));
                        snprintf(buf, 256,"%s",m->mapq_map->description[mapq_idx]);
                        RUN(tld_append(o,buf));
                        RUN(tld_append(o,"',\n"));


                        RUN(tld_append(o,"x: ["));
                        for(int i = 0 ; i < q->len;i++){
                                snprintf(buf, 256,"%d,",i+1);
                                RUN(tld_append(o,buf));
                        }
                        o->len--;
                        RUN(tld_append(o,"],\n"));

                        RUN(tld_append(o,"y: ["));
                        for(int i = 0 ; i < q->len;i++){
                                snprintf(buf, 256,"%f,",mean[i]);
                                RUN(tld_append(o,buf));
                        }
                        o->len--;
                        RUN(tld_append(o,"],\n"));

                        RUN(tld_append(o,"error_y: {\n"));
                        RUN(tld_append(o,"type: 'data',\n"));
                        RUN(tld_append(o,"array : ["));
                        for(int i = 0 ; i < q->len;i++){
                                snprintf(buf, 256,"%f,",stderr[i]);
                                RUN(tld_append(o,buf));
                        }
                        o->len--;
                        RUN(tld_append(o,"],\n"));
                        RUN(tld_append(o,"visible: true\n"));
                        RUN(tld_append(o,"},\n"));

                        RUN(tld_append(o,"type: 'scatter'\n"));
                        RUN(tld_append(o,"}\n"));
                        /* RUN(tld_append(o,"];\n")); */


                        /* RUN(tld_append(o,"y: [")); */
                        /* for(int i = 0 ; i < q->len;i++){ */
                        /*         double mean; */
                        /*         for(int j = 0; j < q->L;j++){ */
                        /*                 snprintf(buf, 256,"%d,",q->data[j][i]); */
                        /*                 RUN(tld_append(o,buf)); */
                        /*         } */
                        /* } */
                        /* o->len--; */
                        /* RUN(tld_append(o,"],\n")); */
                        /* RUN(tld_append(o,"x: Basex,\n")); */
                        /* snprintf(buf, 256,"name: '%s',\n", m->mapq_map->description[mapq_idx]); */
                        /* RUN(tld_append(o,buf)); */
                        /* RUN(tld_append(o,"type: 'box'\n")); */

                        /* RUN(tld_append(o,"};\n")); */
                }
        }

        RUN(tld_append(o,"var qualdata = ["));
        for(int mapq_idx = 0; mapq_idx < 6; mapq_idx++){
                q = m->qual_comp_R1[mapq_idx];
                if(q->n_counts > 0){
                        snprintf(buf, 256,"qualtrace%d,",mapq_idx);
                        RUN(tld_append(o,buf));
                }
        }
        o->len--;
        RUN(tld_append(o,"];\n"));


        RUN(tld_append(o,"var quallayout = {boxmode: 'group'};\n"));

        RUN(tld_append(o,"Plotly.newPlot('qualcomp', qualdata, quallayout);\n"));


/*         var trace1 = { */
/*   y: [0.2, 0.2, 0.6, 1.0, 0.5, 0.4, 0.2, 0.7, 0.9, 0.1, 0.5, 0.3], */
/*   x: x, */
/*   name: 'kale', */
/*   marker: {color: '#3D9970'}, */
/*   type: 'box' */
/* }; */

        /* } */
        /* {} */


        /*         RUN(tld_append(o,"x: [")); */
        /*         snprintf(buf, 256,"%d",1); */
        /*         RUN(tld_append(o,buf)); */
        /*         for(int j = 1 ; j < seq_comp->len;j++){ */
        /*                 snprintf(buf, 256,",%d",j+1); */
        /*                 RUN(tld_append(o,buf)); */
        /*         } */
        /*         RUN(tld_append(o,"],\n")); */
        /* } */
        RUN(tld_append(o,"</script>\n"));

        gfree(mean);
        gfree(stderr);
        return OK;
ERROR:
        gfree(mean);
        gfree(stderr);
        return FAIL;
}

int base_composition_section(tld_strbuf *o, struct metrics *m)
{
        struct seq_composition* seq_comp = NULL;
        char buf[256];
        char nuc[5] = "ACGTN";

        /* if(m->is_aligned == 0){ */

        /* } */

        for(int mapq_idx = 0; mapq_idx < 6; mapq_idx++){
                seq_comp = m->seq_comp_R1[mapq_idx];
                if(seq_comp->n_counts > 0){

                        RUN(tld_append(o, "<h2>Base composition: "));
                        snprintf(buf, 256,"%s",  m->mapq_map->description[mapq_idx] );
                        RUN(tld_append(o, buf));

                        if(seq_comp->n_counts == 1){
                                snprintf(buf, 256," (%d sequence)",  seq_comp->n_counts);
                        }else{
                                snprintf(buf, 256," (%d sequences)",  seq_comp->n_counts);
                        }
                        RUN(tld_append(o, buf));

                        RUN(tld_append(o, "</h2>\n"));

                        RUN(tld_append(o, "<div id=\""));
                        snprintf(buf, 256,"basecomp%d", mapq_idx);
                        RUN(tld_append(o, buf));
                        RUN(tld_append(o, "\" style=\"width:100%;max-width:1400px\"></div>\n"));
                        RUN(tld_append(o,"<script>\n"));

                        for(int i = 0; i < seq_comp->L;i++){
                                snprintf(buf, 256,"var trace%d_%d = {\n",mapq_idx,i );
                                RUN(tld_append(o,buf));
                                RUN(tld_append(o,"x: ["));
                                snprintf(buf, 256,"%d",1);
                                RUN(tld_append(o,buf));
                                for(int j = 1 ; j < seq_comp->len;j++){
                                        snprintf(buf, 256,",%d",j+1);
                                        RUN(tld_append(o,buf));
                                }
                                RUN(tld_append(o,"],\n"));
                                RUN(tld_append(o,"y: ["));
                                snprintf(buf, 256,"%d", seq_comp->data[i][0]);
                                RUN(tld_append(o,buf));
                                for(int j = 1 ; j < seq_comp->len;j++){
                                        snprintf(buf, 256,",%d",seq_comp->data[i][j]);
                                        RUN(tld_append(o,buf));
                                }
                                RUN(tld_append(o,"],\n"));
                                LOG_MSG("L? %d", seq_comp->L);
                                if(seq_comp->L < 20){
                                        snprintf(buf, 256,"name: '%c',\n", nuc[i]);
                                        RUN(tld_append(o,buf));
                                }
                                RUN(tld_append(o,"type: 'bar'\n"));
                                RUN(tld_append(o,"};\n"));
                        }

                        RUN(tld_append(o,"var "));
                        snprintf(buf, 256,"base_comp_data%d", mapq_idx);
                        RUN(tld_append(o,buf ));
                        RUN(tld_append(o," = ["));
                        /* RUN(tld_append(o,"var base_comp_data = [")); */
                        snprintf(buf, 256,"trace%d_%d",mapq_idx,0);
                        RUN(tld_append(o,buf ));
                        for(int i = 1; i < seq_comp->L;i++){
                                snprintf(buf, 256,",trace%d_%d",mapq_idx,i);
                                RUN(tld_append(o,buf ));

                        }
                        RUN(tld_append(o,"];\n"));

                        RUN(tld_append(o,"var base_comp_layout = {barmode: 'stack'};\n"));
                        RUN(tld_append(o,"Plotly.newPlot('"));
                        snprintf(buf, 256,"basecomp%d", mapq_idx);
                        RUN(tld_append(o, buf));
                        RUN(tld_append(o,"', "));

                        snprintf(buf, 256,"base_comp_data%d", mapq_idx);
                        RUN(tld_append(o,buf ));
                        RUN(tld_append(o,", base_comp_layout);\n"));

                        RUN(tld_append(o,"</script>\n"));
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int create_report(struct metrics *m, struct samstat_param *p)
{
        tld_strbuf* out = NULL;
        FILE* f_ptr = NULL;
        RUN(tld_strbuf_alloc(&out, 1024));

        RUN(report_header(out));

        RUN(mapping_quality_overview_section(out,m));
        RUN(base_composition_section(out, m));
        RUN(base_quality_section(out,m));

        RUN(report_footer(out));

        /* if(p->outfile){ */

        f_ptr = fopen("test.html", "w");
        fprintf(f_ptr,"%s", TLD_STR(out));
        fclose(f_ptr);
        /* } */
        tld_strbuf_free(out);
        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        return FAIL;
}
