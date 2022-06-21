#include "core/tld-core.h"
#include "misc/misc.h"
#include "seq/tld-seq.h"
#include "stats/basic.h"
#include "string/str.h"
#include "tld.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#define REPORT_IMPORT
#include "report.h"

/* static int add_count_data(tld_strbuf *o,char *name, char* label,char* color,char* type, uint32_t *x, uint32_t* y,  int len); */



static int print_title(tld_strbuf *o,struct samstat_param *p,int id);

static int mapping_quality_overview_section(tld_strbuf *o, struct metrics *m);
static int base_composition_section(tld_strbuf *o, struct metrics *m, int read);
/* static int base_quality_section(tld_strbuf *o, struct metrics *m); */
static int base_quality_section(tld_strbuf *o, struct metrics *m, int read);

static int error_composition_section(tld_strbuf *o, struct metrics *m, int read);
static int mismatch_composition_section(tld_strbuf *o, struct metrics *m, int read);
static int ins_composition_section(tld_strbuf *o, struct metrics *m, int read);
static int del_composition_section(tld_strbuf *o, struct metrics *m, int read);

static int add_count_data(tld_strbuf *o,char *name, char* label,char* color,char* type, int vis, uint32_t *x, uint32_t* y,  int len);
static int report_header(tld_strbuf *out_buffer);
static int report_footer(tld_strbuf *out_buffer);


int create_report(struct metrics *m, struct samstat_param *p, int id)
{
        tld_strbuf* out = NULL;
        FILE* f_ptr = NULL;

        if(p){
                LOG_MSG("All good");
        }

        RUN(tld_strbuf_alloc(&out, 1024));

        RUN(report_header(out));


        RUN(print_title(out, p,id));

        RUN(mapping_quality_overview_section(out,m));
        if(m->is_paired){
                RUN(base_composition_section(out, m,1));
                RUN(base_composition_section(out, m,2));
        }else{
                RUN(base_composition_section(out, m,0));
        }
        if(m->is_paired){
                RUN(base_quality_section(out,m,1));
                RUN(base_quality_section(out,m,2));
        }else{
                RUN(base_quality_section(out,m,0));
        }
        if(m->is_paired){
                RUN(error_composition_section(out, m,1));
                RUN(error_composition_section(out, m,2));
        }else{
                RUN(error_composition_section(out, m,0));
        }
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

int print_title(tld_strbuf *o,struct samstat_param *p , int id)

{
        char* n = NULL;
        RUN(tlfilename(p->infile[id], &n));
        RUN(tld_append(o,"<h1>"));
        RUN(tld_append(o,n));
        RUN(tld_append(o,"</h1>\n"));

        MFREE(n);
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

        /* RUN(tld_append(o, "<h2>Mapping Stats</h2>\n")); */
        /* RUN(tld_append(o, "<div id=\"MappingStats\" style=\"width:100%;max-width:1400px\"></div>\n")); */

        /* RUN(tld_append(o,"<script>\n")); */
        /* RUN(tld_append(o, "var mapping_data = [{ \n")); */
        /* RUN(tld_append(o, "type: \"pie\",\n")); */
        /* RUN(tld_append(o, "values: [")); */
        /* snprintf(buf, 256, "%d",m->seq_comp_R1[0]->n_counts); */
        /* RUN(tld_append(o,buf)); */
        /* for(int i = 1; i < m->n_mapq_bins;i++){ */
        /*         snprintf(buf, 256, ",%d",m->seq_comp_R1[i]->n_counts); */
        /*         RUN(tld_append(o,buf)); */
        /* } */
        /* RUN(tld_append(o, "],")); */
        /* RUN(tld_append(o, "labels: [")); */
        /* snprintf(buf, 256, "\"%s (%d sequences)\"",m->mapq_map->description[0],m->seq_comp_R1[0]->n_counts); */
        /* RUN(tld_append(o,buf)); */
        /* for(int i = 1; i < m->mapq_map->n_bin;i++){ */
        /*         snprintf(buf, 256, ",\"%s (%d sequences)\"",m->mapq_map->description[i],m->seq_comp_R1[i]->n_counts ); */
        /*         RUN(tld_append(o,buf)); */
        /* } */
        /* RUN(tld_append(o, "],\n")); */

        /* RUN(tld_append(o, "textinfo: \"label+percent\",\n")); */
        /* RUN(tld_append(o, "textposition: \"outside\",\n")); */
        /* RUN(tld_append(o, "automargin: true\n")); */

        /* RUN(tld_append(o, "}]\n")); */

        /* RUN(tld_append(o, "var mapping_layout = {\n")); */
        /* RUN(tld_append(o, "height: 600,\n")); */
        /* RUN(tld_append(o, "width: 600,\n")); */
        /* RUN(tld_append(o, "margin: {\"t\": 0, \"b\": 0, \"l\": 0, \"r\": 0},\n")); */
        /* RUN(tld_append(o, "showlegend: true\n")); */
        /* RUN(tld_append(o, "}\n")); */

        /* RUN(tld_append(o, "Plotly.newPlot('MappingStats', mapping_data, mapping_layout)\n")); */

        /* RUN(tld_append(o,"</script>\n")); */

        /* TABLE  */

        RUN(tld_append(o, "<div id=\"MappingStatsTable\" style=\"width:100%;max-width:700px;height:240px;\"></div>\n"));
        RUN(tld_append(o,"<script>\n"));
        RUN(tld_append(o, "var mappingtableval = [\n"));
        RUN(tld_append(o, "["));
        for(int i = 0; i < m->mapq_map->n_bin;i++){
                snprintf(buf, 256, "\"%s\",",m->mapq_map->description[i]);
                RUN(tld_append(o,buf));
        }
        snprintf(buf, 256, "\"Total\",");
        RUN(tld_append(o,buf));
        o->len--;
        RUN(tld_append(o, "],"));
        uint32_t total = 0;
        RUN(tld_append(o, "["));
        for(int i = 0; i < m->mapq_map->n_bin;i++){
                total += m->seq_comp_R1[i]->n_counts;
                snprintf(buf, 256, "%d,", m->seq_comp_R1[i]->n_counts);
                RUN(tld_append(o,buf));
        }
        snprintf(buf, 256, "%d,", total);
        RUN(tld_append(o,buf));
        o->len--;
        RUN(tld_append(o, "],"));
         RUN(tld_append(o, "["));
        for(int i = 0; i < m->mapq_map->n_bin;i++){
                snprintf(buf, 256, "%.5f,", (double) m->seq_comp_R1[i]->n_counts / (double)total);
                RUN(tld_append(o,buf));
        }
        snprintf(buf, 256, "%.5f,", (double) total / (double )total);
        RUN(tld_append(o,buf));
        o->len--;
        RUN(tld_append(o, "]]\n"));

        RUN(tld_append(o, "var mappingtable_data = [{\n"));
        RUN(tld_append(o, "type: 'table',\n"));
        RUN(tld_append(o, "header: {\n"));
        RUN(tld_append(o, "values: [[\"<b>Mapping Quality Bin</b>\"], [\"<b>Number of Sequences</b>\"], [\"<b>Percentage</b>\"]],\n"));
        RUN(tld_append(o, "align: \"center\",\n"));
        RUN(tld_append(o, "line: {width: 1, color: 'black'},\n"));
        RUN(tld_append(o, "fill: {color: \"grey\"},\n"));
        RUN(tld_append(o, "font: {size: 14, color: \"white\"}\n"));
        RUN(tld_append(o, "},\n"));
        RUN(tld_append(o, "cells: {\n"));
        RUN(tld_append(o, "height: 40,\n"));
        RUN(tld_append(o, "values: mappingtableval,\n"));
        RUN(tld_append(o, "format : [\"\",\",\",\".1%\"],\n"));
        RUN(tld_append(o, "align: [\"center\",\"right\"],\n"));
        RUN(tld_append(o, "line: {color: \"black\", width: 1},\n"));
        RUN(tld_append(o, "font: {size: 14, color: [\"black\"]}\n"));
        RUN(tld_append(o, "}\n"));
        RUN(tld_append(o, "}]\n"));

        RUN(tld_append(o, "var mapping_table_layout = {\n"));
        RUN(tld_append(o, "margin: {\"t\": 0, \"b\": 0, \"l\": 0, \"r\": 0},\n"));
        RUN(tld_append(o, "showlegend: false\n"));
        RUN(tld_append(o, "}\n"));


        RUN(tld_append(o, "Plotly.newPlot('MappingStatsTable', mappingtable_data,mapping_table_layout);\n"));

        RUN(tld_append(o, "</script>\n"));

        return OK;
ERROR:
        return FAIL;
}

int base_quality_section(tld_strbuf *o, struct metrics *m, int read)
{
        struct qual_composition* q = NULL;
        char buf[256];
        char target[16];
        double* mean = NULL;
        double* stderr = NULL;
        uint32_t t = 0;
        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                if(!m->is_paired){
                        /* snprintf(target, 16,"qualcomp"); */
                        /* RUN(tld_append(o, "<h2>Base quality distribution</h2>\n")); */
                        q = m->qual_comp_R1[mapq_idx];
                }else{
                        if(read == 1){
                                /* snprintf(target, 16,"qualcompR1"); */
                                /* RUN(tld_append(o, "<h2>Base quality distribution R1</h2>\n")); */
                                q = m->qual_comp_R1[mapq_idx];
                        }else if(read == 2){
                                /* snprintf(target, 16,"qualcompR2"); */
                                /* RUN(tld_append(o, "<h2>Base quality distribution R2</h2>\n")); */
                                q = m->qual_comp_R2[mapq_idx];
                        }else{
                                ERROR_MSG("Samstat does not support protocols producing more than 2 reads.");
                        }
                }
                t += q->n_counts;
        }
        if(t == 0){
                return OK;
        }
        if(!m->is_paired){
                snprintf(target, 16,"qualcomp");
                RUN(tld_append(o, "<h2>Base quality distribution</h2>\n"));
                q = m->qual_comp_R1[0];
        }else{
                if(read == 1){
                        snprintf(target, 16,"qualcompR1");
                        RUN(tld_append(o, "<h2>Base quality distribution R1</h2>\n"));
                        q = m->qual_comp_R1[0];
                }else if(read == 2){
                        snprintf(target, 16,"qualcompR2");
                        RUN(tld_append(o, "<h2>Base quality distribution R2</h2>\n"));
                        q = m->qual_comp_R2[0];
                }else{
                        ERROR_MSG("Samstat does not support protocols producing more than 2 reads.");
                }
        }
        /* kind of important */
        galloc(&mean, q->len);
        galloc(&stderr, q->len);

        RUN(tld_append(o, "<div id=\""));
        RUN(tld_append(o, target));
        RUN(tld_append(o, "\" style=\"width:100%;max-width:1400px\"></div>\n"));

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

        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                switch (read) {
                case 0:
                case 1:
                        q = m->qual_comp_R1[mapq_idx];
                        break;
                case 2:
                        q = m->qual_comp_R2[mapq_idx];
                        break;
                }

                if(q->n_counts > 0){
                        /* calc mean */
                        for(int i = 0 ; i < q->len;i++){
                                double total = 0.0;

                                double n = 0.0;
                                for(int j = 0; j < q->L;j++){
                                        total += j * q->data[i][j];
                                        n+= q->data[i][j];
                                }
                                if(n == 0){
                                        mean[i] = 0.0;
                                        stderr[i] = 0.0;
                                        /* LOG_MSG("should not happen : %f total: %f", n,total); */
                                }else{
                                        mean[i] = total / n;
                                        total = 0.0;
                                        n = 0.0;
                                        for(int j = 0; j < q->L;j++){
                                                double s = mean[i] - (double) j ;
                                                total += q->data[i][j] *  s * s;
                                                n+= q->data[i][j];
                                        }
                                        stderr[i] = sqrt(total / n) / sqrt(n);
                                }
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
                }
        }

        RUN(tld_append(o,"var qualdata = ["));
        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                switch (read) {
                case 0:
                case 1:
                        q = m->qual_comp_R1[mapq_idx];
                        break;
                case 2:
                        q = m->qual_comp_R2[mapq_idx];
                        break;
                }
                /* q = m->qual_comp_R1[mapq_idx]; */
                if(q->n_counts > 0){
                        snprintf(buf, 256,"qualtrace%d,",mapq_idx);
                        RUN(tld_append(o,buf));
                }
        }
        o->len--;
        RUN(tld_append(o,"];\n"));

        RUN(tld_append(o,"var quallayout = {boxmode: 'group',\n"));
        RUN(tld_append(o,"xaxis: {\n"));
        RUN(tld_append(o,"title: 'Length'\n"));
        RUN(tld_append(o,"},\n"));
        RUN(tld_append(o,"yaxis: {\n"));
        RUN(tld_append(o,"title: 'Base Quality'\n"));
        RUN(tld_append(o,"}\n"));
        RUN(tld_append(o,"};\n"));

        RUN(tld_append(o,"Plotly.newPlot('"));
        RUN(tld_append(o,target)) ;
        RUN(tld_append(o,"', qualdata, quallayout);\n"));

        RUN(tld_append(o,"</script>\n"));

        gfree(mean);
        gfree(stderr);
        return OK;
ERROR:
        gfree(mean);
        gfree(stderr);
        return FAIL;
}

/* helper function to write arrays */

int base_composition_section(tld_strbuf *o, struct metrics *m, int read)
{
        struct seq_composition* seq_comp = NULL;
        char buf[256];
        char name[256];
        char nuc[5] = "ACGTN";
        char target[16];
        /* if(m->is_aligned == 0){ */

        /* } */
        if(!m->is_paired){
                snprintf(target, 16,"basecomp");
                RUN(tld_append(o, "<h2>Base composition</h2>"));

        }else{
                if(read == 1){
                        snprintf(target, 16,"basecompR1");
                        RUN(tld_append(o, "<h2>Base composition R1</h2>"));

                }else if(read == 2){
                        snprintf(target, 16,"basecompR2");
                        RUN(tld_append(o, "<h2>Base composition R2</h2>"));
                }else{
                        ERROR_MSG("Samstat does not support protocols producing more than 2 reads.");
                }
        }


        RUN(tld_append(o, "<div id=\""));
        /* snprintf(buf, 256,"basecomp"); */
        RUN(tld_append(o, target));
        RUN(tld_append(o, "\" style=\"width:100%;max-width:1400px\"></div>\n"));
        RUN(tld_append(o,"<script>\n"));
        int vis =1;

        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                switch (read) {
                case 0:
                case 1:
                        seq_comp = m->seq_comp_R1[mapq_idx];
                        break;
                case 2:
                        seq_comp = m->seq_comp_R2[mapq_idx];
                        break;
                }

                if(seq_comp->n_counts > 0){

                        for(int i = 0; i < seq_comp->L;i++){
                                snprintf(buf, 256,"trace%d_%d",mapq_idx,i );
                                snprintf(name, 256,"%c",nuc[i]);
                                RUN(add_count_data(
                                            o,
                                            buf,
                                            name,
                                            NULL,
                                            "bar",
                                            vis,
                                            NULL,
                                            seq_comp->data[i],
                                            seq_comp->len
                                            ));
                        }
                        vis = 0;
                }
        }
        RUN(tld_append(o,"var "));
        snprintf(buf, 256,"base_comp_data");
        RUN(tld_append(o,buf ));
        RUN(tld_append(o," = ["));
        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                switch (read) {
                case 0:
                case 1:
                        seq_comp = m->seq_comp_R1[mapq_idx];
                        break;
                case 2:
                        seq_comp = m->seq_comp_R2[mapq_idx];
                        break;
                }


                if(seq_comp->n_counts > 0){
                        for(int i = 0; i < seq_comp->L;i++){
                                snprintf(buf, 256,"trace%d_%d,",mapq_idx,i);
                                RUN(tld_append(o,buf ));
                        }
                }

        }
        o->len--;
        RUN(tld_append(o,"];\n"));
        RUN(tld_append(o,"var base_comp_layout = {barmode: 'stack'};\n"));
        RUN(tld_append(o,"Plotly.newPlot('"));
        /* snprintf(buf, 256,"basecomp"); */
        RUN(tld_append(o, target));
        RUN(tld_append(o,"', "));

        snprintf(buf, 256,"base_comp_data,");
        RUN(tld_append(o,buf ));
        /* RUN(tld_append(o,"base_comp_layout,")); */

        RUN(tld_append(o,"{barmode: 'stack',\n"));
        RUN(tld_append(o,"title: 'Base Composition',\n"));

        RUN(tld_append(o,"xaxis: {\n"));
        RUN(tld_append(o,"title: 'Length'\n"));
        RUN(tld_append(o,"},\n"));
        RUN(tld_append(o,"yaxis: {\n"));
        RUN(tld_append(o,"title: 'Counts'\n"));
        RUN(tld_append(o,"},\n"));

        RUN(tld_append(o,"updatemenus: [{\n"));
        RUN(tld_append(o,"x: 0.5,\n"));
        RUN(tld_append(o,"xanchor: 'center',\n"));

        RUN(tld_append(o,"y: 1.09,\n"));
        RUN(tld_append(o,"yanchor: 'top',\n"));
        /* RUN(tld_append(o,"y: 1,\n")); */
        /* RUN(tld_append(o,"yanchor: 'top',\n")); */
        RUN(tld_append(o,"buttons: [\n"));
        for(int vis = 0; vis < m->n_mapq_bins; vis++){
                seq_comp = m->seq_comp_R1[vis];
                if(seq_comp->n_counts > 0){
                        RUN(tld_append(o,"{\n"));
                        RUN(tld_append(o,"method: 'restyle',\n"));
                        RUN(tld_append(o,"args: ['visible', [\n"));
                        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                                seq_comp = m->seq_comp_R1[mapq_idx];
                                if(seq_comp->n_counts > 0){
                                        for(int i = 0; i < seq_comp->L;i++){
                                                if(vis == mapq_idx){
                                                        RUN(tld_append(o,"true,"));
                                                }else{
                                                        RUN(tld_append(o,"false,"));
                                                }
                                        }
                                }

                        }
                        o->len--;
                        RUN(tld_append(o,"]],\n"));
                        snprintf(buf,256,"label: '%s'\n",m->mapq_map->description[vis]);
                        RUN(tld_append(o,buf));
                        RUN(tld_append(o,"},\n"));
                }
        }
        o->len--;
        RUN(tld_append(o,"]}]}\n"));

        RUN(tld_append(o,");\n"));
        RUN(tld_append(o,"</script>\n"));
        return OK;
ERROR:
        return FAIL;
}

int error_composition_section(tld_strbuf *o, struct metrics *m, int read)
{
        if(m->is_aligned){

                if(!m->is_paired){
                        RUN(tld_append(o, "<h2>Error composition</h2>"));

                }else{
                        if(read == 1){
                          RUN(tld_append(o, "<h2>Error composition R1</h2>"));
                        }else if(read == 2){
                                RUN(tld_append(o, "<h2>Error composition R2</h2>"));
                        }else{
                                ERROR_MSG("Samstat does not support protocols producing more than 2 reads.");
                        }
                }

                        /* RUN(tld_append(o, "<h2>Error composition: </h2>")); */
                RUN(mismatch_composition_section(o,m,read));
                RUN(ins_composition_section(o,m,read));
                RUN(del_composition_section(o,m,read));
        }
        return OK;
ERROR:
        return FAIL;
}

int mismatch_composition_section(tld_strbuf *o, struct metrics *m, int read)
{
        struct error_composition* e = NULL;
        char buf[256];
        char name[256];
        char target[16];
        char nuc[5] = "ACGTN";

        if(!m->is_paired){
                snprintf(target, 16,"miserrcomp");
        }else{
                if(read == 1){
                        snprintf(target, 16,"miserrcompR1");
                }else if(read == 2){
                        snprintf(target, 16,"miserrcompR2");
                }else{
                        ERROR_MSG("Samstat does not support protocols producing more than 2 reads.");
                }
        }

        RUN(tld_append(o, "<div id=\""));
        /* snprintf(buf, 256,"miserrcomp"); */
        RUN(tld_append(o, target));
        RUN(tld_append(o, "\" style=\"width:100%;max-width:1400px\"></div>\n"));
        RUN(tld_append(o,"<script>\n"));
        int vis = 1;
        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                switch (read) {
                case 0:
                case 1:
                        /* seq_comp = m->seq_comp_R1[mapq_idx]; */
                        e = m->error_comp_R1[mapq_idx];
                        break;
                case 2:
                        e = m->error_comp_R2[mapq_idx];
                        /* seq_comp = m->seq_comp_R2[mapq_idx]; */
                        break;
                }

                /* e = m->error_comp_R1[mapq_idx]; */
                /* LOG_MSG("%s %d",m->mapq_map->description[mapq_idx],e->n_mis + e->n_del + e->n_ins ); */
                if(e->n_mis + e->n_del + e->n_ins  > 0){

                        /* double  total = (double) e->n_mis ;//+ e->n_del + e->n_ins; */
                        /* double  total = (double) m->seq_comp_R1[mapq_idx]->n_counts; */
                        /* uint32_t sanity = 0; */
                        /* for(int i = 0; i < e->L;i++){ */
                        /*         for(int j = 0; j < e->len;j++){ */
                                        /* e->mis[i][j] = round( (double) e->mis[i][j] * 1000000.0 / total); */
                                        /* sanity+= e->mis[i][j]; */
                        /*         } */
                        /* } */
                        /* LOG_MSG("Sanity: %d total %d", sanity, total); */
                        for(int i = 0; i < e->L;i++){
                                snprintf(buf, 256,"errtrace%d_%d",mapq_idx,i );
                                snprintf(name, 256,"%c",nuc[i]);
                                RUN(add_count_data(
                                            o,
                                            buf,
                                            name,
                                            NULL,
                                            "bar",
                                            vis,
                                            NULL,
                                            e->mis[i],
                                            e->len
                                            ));
                        }
                        vis = 0;
                }
        }
        RUN(tld_append(o,"var "));
        snprintf(buf, 256,"err_comp_data");
        RUN(tld_append(o,buf ));
        RUN(tld_append(o," = ["));
        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                switch (read) {
                case 0:
                case 1:
                        /* seq_comp = m->seq_comp_R1[mapq_idx]; */
                        e = m->error_comp_R1[mapq_idx];
                        break;
                case 2:
                        e = m->error_comp_R2[mapq_idx];
                        /* seq_comp = m->seq_comp_R2[mapq_idx]; */
                        break;
                }

                /* e = m->error_comp_R1[mapq_idx]; */
                /* LOG_MSG("%s %d",m->mapq_map->description[mapq_idx],e->n_mis + e->n_del + e->n_ins ); */
                if(e->n_mis + e->n_del + e->n_ins  > 0){


                        for(int i = 0; i < e->L;i++){
                                snprintf(buf, 256,"errtrace%d_%d,",mapq_idx,i );
                                RUN(tld_append(o,buf ));
                        }
                }
        }
        o->len--;


        RUN(tld_append(o,"];\n"));

        RUN(tld_append(o,"Plotly.newPlot('"));
        /* snprintf(buf, 256,"miserrcomp"); */
        RUN(tld_append(o, target));
        RUN(tld_append(o,"', "));

        snprintf(buf, 256,"err_comp_data,");
        RUN(tld_append(o,buf ));
        /* RUN(tld_append(o,"base_comp_layout,")); */

        RUN(tld_append(o,"{barmode: 'stack',\n"));
        RUN(tld_append(o,"title: 'Distribution of Mismatches',\n"));

        RUN(tld_append(o,"xaxis: {\n"));
        RUN(tld_append(o,"title: 'Length'\n"));
        RUN(tld_append(o,"},\n"));
        RUN(tld_append(o,"yaxis: {\n"));
        RUN(tld_append(o,"title: 'Counts'\n"));
        RUN(tld_append(o,"},\n"));


        RUN(tld_append(o,"updatemenus: [{\n"));
        RUN(tld_append(o,"x: 0.5,\n"));
        RUN(tld_append(o,"xanchor: 'center',\n"));

        RUN(tld_append(o,"y: 1.09,\n"));
        RUN(tld_append(o,"yanchor: 'top',\n"));
        RUN(tld_append(o,"buttons: [\n"));
        for(int vis = 0; vis < m->n_mapq_bins; vis++){
                e = m->error_comp_R1[vis];
                /* LOG_MSG("%s %d",m->mapq_map->description[mapq_idx],e->n_mis + e->n_del + e->n_ins ); */
                if(e->n_mis + e->n_del + e->n_ins  > 0){
                        RUN(tld_append(o,"{\n"));
                        RUN(tld_append(o,"method: 'restyle',\n"));
                        RUN(tld_append(o,"args: ['visible', [\n"));
                        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                                e = m->error_comp_R1[mapq_idx];
                                if(e->n_mis + e->n_del + e->n_ins  > 0){
                                        for(int i = 0; i < e->L;i++){
                                                if(vis == mapq_idx){
                                                        RUN(tld_append(o,"true,"));
                                                }else{
                                                        RUN(tld_append(o,"false,"));
                                                }
                                        }
                                }
                        }
                        o->len--;
                        RUN(tld_append(o,"]],\n"));
                        snprintf(buf,256,"label: '%s'\n",m->mapq_map->description[vis]);
                        RUN(tld_append(o,buf));
                        RUN(tld_append(o,"},\n"));
                }
        }
        o->len--;
        RUN(tld_append(o,"]}]}\n"));

        RUN(tld_append(o,");\n"));
        RUN(tld_append(o,"</script>\n"));
        return OK;
ERROR:
        return FAIL;
}

int ins_composition_section(tld_strbuf *o, struct metrics *m, int read)
{
        struct error_composition* e = NULL;
        char buf[256];
        char name[256];
        char target[16];

        if(!m->is_paired){
                snprintf(target, 16,"inserrcomp");
        }else{
                if(read == 1){
                        snprintf(target, 16,"inserrcompR1");
                }else if(read == 2){
                        snprintf(target, 16,"inserrcompR2");
                }else{
                        ERROR_MSG("Samstat does not support protocols producing more than 2 reads.");
                }
        }

        RUN(tld_append(o, "<div id=\""));
        /* snprintf(buf, 256,"inserrcomp"); */
        RUN(tld_append(o, target));
        RUN(tld_append(o, "\" style=\"width:100%;max-width:1400px\"></div>\n"));
        RUN(tld_append(o,"<script>\n"));

        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                switch (read) {
                case 0:
                case 1:
                        /* seq_comp = m->seq_comp_R1[mapq_idx]; */
                        e = m->error_comp_R1[mapq_idx];
                        break;
                case 2:
                        e = m->error_comp_R2[mapq_idx];
                        /* seq_comp = m->seq_comp_R2[mapq_idx]; */
                        break;
                }


                /* LOG_MSG("%s %d",m->mapq_map->description[mapq_idx],e->n_mis + e->n_del + e->n_ins ); */
                if(e->n_mis + e->n_del + e->n_ins  > 0){
                        snprintf(buf, 256,"inserrtrace%d",mapq_idx);
                        snprintf(name, 256,"%s", m->mapq_map->description[mapq_idx]);
                        RUN(add_count_data(
                                    o,
                                    buf,
                                    name,
                                    NULL,
                                    "scatter",
                                    1,
                                    NULL,
                                    e->ins,
                                    e->len
                                    ));
                }
        }
        RUN(tld_append(o,"var "));
        snprintf(buf, 256,"inserr_comp_data");
        RUN(tld_append(o,buf ));
        RUN(tld_append(o," = ["));
        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                switch (read) {
                case 0:
                case 1:
                        /* seq_comp = m->seq_comp_R1[mapq_idx]; */
                        e = m->error_comp_R1[mapq_idx];
                        break;
                case 2:
                        e = m->error_comp_R2[mapq_idx];
                        /* seq_comp = m->seq_comp_R2[mapq_idx]; */
                        break;
                }

                /* e = m->error_comp_R1[mapq_idx]; */
                /* LOG_MSG("%s %d",m->mapq_map->description[mapq_idx],e->n_mis + e->n_del + e->n_ins ); */
                if(e->n_mis + e->n_del + e->n_ins  > 0){
                        snprintf(buf, 256,"inserrtrace%d,",mapq_idx);
                        RUN(tld_append(o,buf ));
                }
        }
        o->len--;
        RUN(tld_append(o,"];\n"));
        RUN(tld_append(o,"Plotly.newPlot('"));
        /* snprintf(buf, 256,"inserrcomp"); */
        RUN(tld_append(o, target));
        RUN(tld_append(o,"', "));

        snprintf(buf, 256,"inserr_comp_data,");
        RUN(tld_append(o,buf ));
        /* RUN(tld_append(o,"base_comp_layout,")); */

        RUN(tld_append(o,"{boxmode: 'group',\n"));
        RUN(tld_append(o,"title: 'Distribution of Insertions',\n"));

        RUN(tld_append(o,"xaxis: {\n"));
        RUN(tld_append(o,"title: 'Length'\n"));
        RUN(tld_append(o,"},\n"));
        RUN(tld_append(o,"yaxis: {\n"));
        RUN(tld_append(o,"title: 'Counts'\n"));
        RUN(tld_append(o,"}}\n"));
        RUN(tld_append(o,");\n"));
        RUN(tld_append(o,"</script>\n"));
        return OK;
ERROR:
        return FAIL;
}

int del_composition_section(tld_strbuf *o, struct metrics *m, int read)
{
        struct error_composition* e = NULL;
        char buf[256];
        char name[256];
                char target[16];

        if(!m->is_paired){
                snprintf(target, 16,"delerrcomp");
        }else{
                if(read == 1){
                        snprintf(target, 16,"delerrcompR1");
                }else if(read == 2){
                        snprintf(target, 16,"delerrcompR2");
                }else{
                        ERROR_MSG("Samstat does not support protocols producing more than 2 reads.");
                }
        }

        RUN(tld_append(o, "<div id=\""));
        /* snprintf(buf, 256,"delerrcomp"); */
        RUN(tld_append(o, target));
        RUN(tld_append(o, "\" style=\"width:100%;max-width:1400px\"></div>\n"));
        RUN(tld_append(o,"<script>\n"));

        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                switch (read) {
                case 0:
                case 1:
                        /* seq_comp = m->seq_comp_R1[mapq_idx]; */
                        e = m->error_comp_R1[mapq_idx];
                        break;
                case 2:
                        e = m->error_comp_R2[mapq_idx];
                        /* seq_comp = m->seq_comp_R2[mapq_idx]; */
                        break;
                }

                /* e = m->error_comp_R1[mapq_idx]; */
                /* LOG_MSG("%s %d",m->mapq_map->description[mapq_idx],e->n_mis + e->n_del + e->n_ins ); */
                if(e->n_mis + e->n_del + e->n_ins  > 0){
                        snprintf(buf, 256,"delerrtrace%d",mapq_idx);
                        snprintf(name, 256,"%s", m->mapq_map->description[mapq_idx]);
                        RUN(add_count_data(
                                    o,
                                    buf,
                                    name,
                                    NULL,
                                    "scatter",
                                    1,
                                    NULL,
                                    e->del,
                                    e->len
                                    ));
                }
        }
        RUN(tld_append(o,"var "));
        snprintf(buf, 256,"delerr_comp_data");
        RUN(tld_append(o,buf ));
        RUN(tld_append(o," = ["));
        for(int mapq_idx = 0; mapq_idx < m->n_mapq_bins; mapq_idx++){
                switch (read) {
                case 0:
                case 1:
                        /* seq_comp = m->seq_comp_R1[mapq_idx]; */
                        e = m->error_comp_R1[mapq_idx];
                        break;
                case 2:
                        e = m->error_comp_R2[mapq_idx];
                        /* seq_comp = m->seq_comp_R2[mapq_idx]; */
                        break;
                }

                                /* e = m->error_comp_R1[mapq_idx]; */
                /* LOG_MSG("%s %d",m->mapq_map->description[mapq_idx],e->n_mis + e->n_del + e->n_ins ); */
                if(e->n_mis + e->n_del + e->n_ins  > 0){
                        snprintf(buf, 256,"delerrtrace%d,",mapq_idx);
                        RUN(tld_append(o,buf ));
                }
        }
        o->len--;
        RUN(tld_append(o,"];\n"));
        RUN(tld_append(o,"Plotly.newPlot('"));
        /* snprintf(buf, 256,"delerrcomp"); */
        RUN(tld_append(o, target));
        RUN(tld_append(o,"', "));

        snprintf(buf, 256,"delerr_comp_data,");
        RUN(tld_append(o,buf ));
        /* RUN(tld_append(o,"base_comp_layout,")); */

        RUN(tld_append(o,"{boxmode: 'group',\n"));
        RUN(tld_append(o,"title: 'Distribution of Deletions',\n"));

        RUN(tld_append(o,"xaxis: {\n"));
        RUN(tld_append(o,"title: 'Length'\n"));
        RUN(tld_append(o,"},\n"));
        RUN(tld_append(o,"yaxis: {\n"));
        RUN(tld_append(o,"title: 'Counts'\n"));
        RUN(tld_append(o,"}}\n"));
        RUN(tld_append(o,");\n"));
        RUN(tld_append(o,"</script>\n"));
        return OK;
ERROR:
        return FAIL;
}


int add_count_data(tld_strbuf *o,char *name, char* label,char* color,char* type, int vis, uint32_t *x, uint32_t* y,  int len)
{
        ASSERT(y != NULL, "no y");
        char buf[256];
        snprintf(buf, 256,"var %s  = {\n",name);
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

        }
        RUN(tld_append(o,"y: ["));
        snprintf(buf, 256,"%d", y[0]);
        RUN(tld_append(o,buf));
        for(int i = 1 ; i < len;i++){
                snprintf(buf, 256,",%d",y[i]);
                RUN(tld_append(o,buf));
        }
        RUN(tld_append(o,"],\n"));
        if(label){
                snprintf(buf, 256,"name: '%s',\n", label);
                RUN(tld_append(o,buf));
        }
        if(color){
                /* snprintf(buf, 256,"name: %s,\n", label); */
                /* RUN(tld_append(o,buf)); */
        }
        if(type){
                snprintf(buf, 256,"type: '%s',\n", type);
                RUN(tld_append(o,buf));
        }else{
                RUN(tld_append(o,"type: 'scatter',\n"));
        }
        if(!vis){
                RUN(tld_append(o,"visible: false\n"));
        }else{
                RUN(tld_append(o,"visible: true\n"));
        }

        RUN(tld_append(o,"};\n"));
        return OK;
ERROR:
        return FAIL;
}


int report_header(tld_strbuf *out_buffer)
{
        RUN(tld_append(out_buffer, "<!DOCTYPE html>\n"));
        RUN(tld_append(out_buffer, "<html>\n"));
        RUN(tld_append(out_buffer, "<script src=\"https://d3js.org/d3.v4.js\"></script>\n"));
        RUN(tld_append(out_buffer, "<script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>\n"));
        /* RUN(tld_append(out_buffer, "<title>%s</title>\n", "ARFFF"); */

        RUN(tld_append(out_buffer, "<style>\n"));
        //RUN(tld_append(out_buffer, "canvas{\n"));
        //RUN(tld_append(out_buffer, "}\n"));

        RUN(tld_append(out_buffer, "body {\n"));
        RUN(tld_append(out_buffer, "	margin: 0 auto;\n"));
        RUN(tld_append(out_buffer, "	padding: 22px 0;\n"));
        RUN(tld_append(out_buffer, "	font-family: 'Source Sans Pro', sans-serif;\n"));
        RUN(tld_append(out_buffer, "	font-weight: 200;\n"));
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
        RUN(tld_append(out_buffer, "	line-height:140%%;\n"));
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

int report_footer(tld_strbuf *out_buffer)
{
        RUN(tld_append(out_buffer, "<footer>\n"));

        RUN(tld_append(out_buffer, "<div>\n"));
        RUN(tld_append(out_buffer, "<section id=\"about\">\n"));
        RUN(tld_append(out_buffer, "<h3>About</h3>\n"));
        RUN(tld_append(out_buffer, "<p>SAMStat is a basic quality control program for NGS data. It displays information about mapped and unmapped reads to diagnose mapping and other problems.</p> <p> The plots on this page were draw using <a href=\"https://plotly.com/graphing-libraries/\">Plotly</a>.</p>\n"));
        RUN(tld_append(out_buffer, "</section>\n"));
        RUN(tld_append(out_buffer, "<section id=\"contact\">\n"));

        RUN(tld_append(out_buffer, "<h3>Contact</h3>\n"));

        RUN(tld_append(out_buffer, "<p>SAMStat was developed by Timo Lassmann (timolassmann at gmail dot com).</p>\n"));
        RUN(tld_append(out_buffer, "<h3>Please cite:</h3>\n"));

        RUN(tld_append(out_buffer, "<p>Lassmann et al. (2011) \"SAMStat: monitoring biases in next generation sequencing data.\" Bioinformatics doi:10.1093/bioinformatics/btq614 [<a href =\"http://www.ncbi.nlm.nih.gov/pubmed/21088025/\">PMID: 21088025</a>] </p>\n"));
        RUN(tld_append(out_buffer, "</section>\n"));
        RUN(tld_append(out_buffer, "<section id=\"blogroll\">\n"));
        RUN(tld_append(out_buffer, "<h3>Links</h3>\n"));
        RUN(tld_append(out_buffer, "<ul>\n"));
        /* RUN(tld_append(out_buffer, "<li><a href=\"http://telethonkids.org.au/\">Telethon Kids Institute</a></li>\n")); */
        RUN(tld_append(out_buffer, "<li><a href=\"https://samtools.github.io/hts-specs/\">SAM/BAM specifications</a></li>\n"));
        RUN(tld_append(out_buffer, "<li><a href=\"https://github.com/TimoLassmann/samstat/\">samstat</a></li>\n"));
        /* RUN(tld_append(out_buffer, "<li><a href=\"http://code.google.com/p/bedtools/\">BEDtools</a></li>\n")); */
        RUN(tld_append(out_buffer, "<li><a href=\"https://github.com/TimoLassmann/kalign/\">Kalign</a></li>\n"));
        /* RUN(tld_append(out_buffer, "<li><a href=\"http://ngsview.sourceforge.net/\">NGSview</a></li>\n")); */


        RUN(tld_append(out_buffer, "</ul>\n"));
        RUN(tld_append(out_buffer, "</section>\n"));

        RUN(tld_append(out_buffer, "</div>\n"));
        RUN(tld_append(out_buffer, "</footer>\n"));

        RUN(tld_append(out_buffer, "</body>\n"));
        RUN(tld_append(out_buffer, "</html>\n"));
        return OK;
ERROR:
        return FAIL;
}
