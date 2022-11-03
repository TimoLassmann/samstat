#include "tld.h"

#include <sys/stat.h>
#include <time.h>
#include "lst.h"
#include "plot.h"
#include "../param/param.h"
#include "../collect/collect.h"
#include "../config.h"
#define  STAT_REPORT_IMPORT
#include "stat_report.h"

static int samstat_title(tld_strbuf *o, struct samstat_param *p, int id);
static int samstat_file_stat_section(tld_strbuf *o, char *filename);
static int samstat_version_section(tld_strbuf *o);
static int samstat_partial_warning_section(tld_strbuf *o, struct stat_collection *s);
static int samstat_mapping_quality_overview_section(tld_strbuf *o, struct stat_collection *s);
static int samstat_footer(tld_strbuf *out_buffer);

static int samstat_add_large_number(tld_strbuf *o, uint64_t n);

static int write_buffer(tld_strbuf *out_buffer, struct samstat_param *p, int id);

int stat_report(struct stat_collection* s, struct samstat_param *p, int id)
{

        char* n = NULL;
        tld_strbuf* out = NULL;
        RUN(tlfilename(p->infile[id], &n));

        RUN(stat_collection_finalise(s));

        RUN(tld_strbuf_alloc(&out, 1024));
        html_header(out,n);

        RUN(samstat_title(out, p, id));

        RUN(samstat_file_stat_section(out, p->infile[id]));
        RUN(samstat_version_section(out));

        /* LOG_MSG("add mapping stat"); */
        RUN(samstat_mapping_quality_overview_section(out, s));

        if(s->is_partial_report){
                RUN(samstat_partial_warning_section(out, s));
        }

        /* for(int i = 0 ; i < 2;i++){ */
        /*         for(int j = 0; j < s->mapq_map->n_bin;j++){ */
        /*                 fprintf(stderr,"%"PRId64"\t", s->basic_nums[i][j]); */
        /*         } */
        /*         fprintf(stderr," - READ %d\n",i+1); */
        /* } */

        /* LENGTH Distribution  */
        if(!s->n_paired){
                RUN(tld_append(out, "<h2>Length distribution</h2>\n"));
        }else{
                RUN(tld_append(out, "<h2>Length distribution Read 1</h2>\n"));
        }
        plot_add(out,s->len_dist_R1);
        if(s->is_aligned){
                RUN(tld_append(out,"<p>Density plot of read lengths split by mapping quality. If aligned reads are soft clipped, additional length distribution of the aligned reads are shown (indicated by the 'mapped' keyword). </p>"));
        }else{
                RUN(tld_append(out,"<p>Density plot of read lengths. </p>"));
        }
        if(s->n_paired){
                RUN(tld_append(out, "<h2>Length distribution Read 2</h2>\n"));
                plot_add(out,s->len_dist_R2);
                if(s->is_aligned){
                        RUN(tld_append(out,"<p>Density plot of read lengths split by mapping quality. If aligned reads are soft clipped, additional length distribution of the aligned reads are shown (indicated by the 'mapped' keyword). </p>"));
                }else{
                        RUN(tld_append(out,"<p>Density plot of read lengths. </p>"));
                }
        }


        /* Base composition   */
        if(!s->n_paired){
                RUN(tld_append(out, "<h2>Base composition</h2>"));
        }else{
                RUN(tld_append(out, "<h2>Base composition Read 1 </h2>"));
        }

        plot_add(out,s->base_comp_R1);
        RUN(tld_append(out,"<p>Base composition. Use the toggle box to view composition of reads mapped at different mapping qualities. Soft clipped sequences are excluded.</p>"));

        if(s->n_paired){
                RUN(tld_append(out, "<h2>Base composition Read 2 </h2>"));
                plot_add(out,s->base_comp_R2);
                RUN(tld_append(out,"<p>Base composition. Use the toggle box to view composition of reads mapped at different mapping qualities. Soft clipped sequences are excluded.</p>"));

        }


        /* Base quality  */

        if(!s->n_paired){
                RUN(tld_append(out, "<h2>Base quality distribution</h2>\n"));
        }else{
                RUN(tld_append(out, "<h2>Base quality distribution Read 1</h2>\n"));
        }

        plot_add(out,s->qual_comp_R1);
        if(s->is_aligned){
                RUN(tld_append(out,"<p>Base quality distribution split up by mapping quality.</p>"));
        }else{
                RUN(tld_append(out,"<p>Base quality distribution.</p>"));
        }
        if(s->n_paired){
                RUN(tld_append(out, "<h2>Base quality distribution Read 2</h2>\n"));
                plot_add(out,s->qual_comp_R2);
                if(s->is_aligned){
                        RUN(tld_append(out,"<p>Base quality distribution split up by mapping quality.</p>"));
                }else{
                        RUN(tld_append(out,"<p>Base quality distribution.</p>"));
                }
        }



        /* end report if  */
        if(s->collect_end){
                if(!s->n_paired){
                        RUN(tld_append(out, "<h2>Base composition at end of reads</h2>"));
                }else{
                        RUN(tld_append(out, "<h2>Base composition at end of Read 1 </h2>"));
                }

                plot_add(out,s->base_comp_R1_end);
                RUN(tld_append(out,"<p>Base compostions at the end of the reads (0), the second last base (1) and so forth. Use the toggle box to view composition of reads mapped at different mapping qualities. Soft clipped sequences are excluded.</p>"));

                if(s->n_paired){
                        RUN(tld_append(out, "<h2>Base composition at end of Read 2 </h2>"));
                        plot_add(out,s->base_comp_R2_end);
                        RUN(tld_append(out,"<p>Base compostions at the end of the reads (0), the second last base (1) and so forth. Use the toggle box to view composition of reads mapped at different mapping qualities. Soft clipped sequences are excluded.</p>"));

                }

                if(!s->n_paired){
                        RUN(tld_append(out, "<h2>Base quality distribution at end of reads</h2>\n"));
                }else{
                        RUN(tld_append(out, "<h2>Base quality distribution at end of Read 1</h2>\n"));
                }
                /* RUN(tld_append(out, "<h2>Base quality distribution at end of Read 1</h2>\n")); */
                plot_add(out,s->qual_comp_R1_end);

                RUN(tld_append(out,"<p>Base quality distribution at the last (0), second last (1) and so forth base.</p>"));

                if(s->n_paired){
                        RUN(tld_append(out, "<h2>Base quality distribution at end of Read 2</h2>\n"));
                        plot_add(out,s->qual_comp_R2_end);
                        RUN(tld_append(out,"<p>Base quality distribution at the last (0), second last (1) and so forth base.</p>"));
                }
        }

        /* norm  */
        if(s->norm_len_plot){
                if(!s->n_paired){
                        RUN(tld_append(out, "<h2>Base composition along the read</h2>"));
                }else{
                        RUN(tld_append(out, "<h2>Base composition along read 1</h2>"));
                }

                plot_add(out,s->base_comp_len_norm_R1);
                RUN(tld_append(out,"<p>Base compostion; reads are split into bins each representing 1% of the read length.</p>"));
                if(s->n_paired){
                        RUN(tld_append(out, "<h2>Base composition along read 2 </h2>"));

                        plot_add(out,s->base_comp_len_norm_R2);
                        RUN(tld_append(out,"<p>Base compostion; reads are split into bins each representing 1% of the read length.</p>"));

                }


                if(!s->n_paired){
                        RUN(tld_append(out, "<h2>Base quality distibution along the read</h2>"));
                }else{
                        RUN(tld_append(out, "<h2>Base quality distibution along read 1</h2>"));
                }
                plot_add(out,s->qual_comp_len_norm_R1);
                RUN(tld_append(out,"<p>Base quality distibution; reads are split into bins each representing 1% of the read length.</p>"));
                if(s->n_paired){
                        RUN(tld_append(out, "<h2>Base quality distibution along read 2 </h2>"));


                        plot_add(out,s->qual_comp_len_norm_R2);
                        RUN(tld_append(out,"<p>Base quality distibution; reads are split into bins each representing 1% of the read length.</p>"));
                }
        }





        if(s->is_aligned){
                if(!s->n_paired){
                        RUN(tld_append(out, "<h2>Error composition</h2>"));
                }else{
                        RUN(tld_append(out, "<h2>Error composition Read 1</h2>"));
                }
                plot_add(out,s->mis_R1);
                plot_add(out,s->ins_R1);
                plot_add(out,s->del_R1);
                RUN(tld_append(out,"<p>Distibution of mismatches, insertions and deletions along the read length. Note: if no 'MD' tag is present in the SAM/BAM file no mismatch profile is displayed.</p>"));

                if(s->n_paired){
                        RUN(tld_append(out, "<h2>Error composition Read 2</h2>"));
                        plot_add(out,s->mis_R2);

                        plot_add(out,s->ins_R2);

                        plot_add(out,s->del_R2);

                        RUN(tld_append(out,"<p>Distibution of mismatches, insertions and deletions along the read length. Note: if no 'MD' tag is present in the SAM/BAM file no mismatch profile is displayed.</p>"));

                }
        }

        RUN(samstat_footer(out));
        RUN(html_end(out));

        RUN(write_buffer(out, p,id));

        tld_strbuf_free(out);
        MFREE(n);
        return OK;
ERROR:
        return FAIL;
}

int write_buffer(tld_strbuf *out_buffer, struct samstat_param *p, int id)
{
        char buffer[1024];
        FILE* f_ptr = NULL;
        if(p->outdir){
                char* filename = NULL;
                tlfilename(p->infile[id], &filename);
                snprintf(buffer,1024,"%s/%s.samstat.html",p->outdir, filename);
                MFREE(filename);
        }else{
                snprintf(buffer,1024,"%s.samstat.html", p->infile[id]);
        }
        LOG_MSG("Writing to: %s", buffer);
        f_ptr = fopen(buffer, "w");
        ASSERT(f_ptr != NULL,"Failed to open output file");

        fprintf(f_ptr,"%s", TLD_STR(out_buffer));
        fclose(f_ptr);
        return OK;
ERROR:
        return FAIL;
}


int samstat_title(tld_strbuf *o,struct samstat_param *p , int id)
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

int samstat_file_stat_section(tld_strbuf *o, char* filename)
{
        struct stat buf;
        int local_ret;
        char buffer[256];
        struct tm newtime;

        local_ret= stat ( filename, &buf );
        if ( local_ret == 0 ){

                localtime_r(&buf.st_mtime, &newtime);
                tld_append(o, "<p>File size: ");
                snprintf(buffer, 256,"%"PRId64"", (uint64_t)buf.st_size);
                tld_append(o, buffer);
                tld_append(o, " bytes, created ");
                if(!strftime(buffer, 256, "%F %H:%M:%S", &newtime )){
                        ERROR_MSG("Could not make time string.");
                }

                tld_append(o, buffer);
                tld_append(o, "</p>");
                /* strftime(time_string, 200, "%F %H:%M:%S\t", ptr); */
                /* sprintf(buffer,"size:%lld bytes, created %s",(long long)buf.st_size, time_string); */
        }else{
                ERROR_MSG("Failed getting stats for file:%s\n",filename );
        }


        return OK;
ERROR:
        return FAIL;
}

int samstat_version_section(tld_strbuf *o)
{
        char buffer[256];
        RUN(tld_append(o,"<p>(report produced by: "));
        snprintf(buffer, 256,"%s version %s",PACKAGE_NAME,PACKAGE_VERSION);
        RUN(tld_append(o,buffer));
        RUN(tld_append(o,")</p>"));
        return OK;
ERROR:
        return FAIL;
}

int samstat_partial_warning_section(tld_strbuf *o , struct stat_collection*s)
{
        char buf[256];
        snprintf(buf, 256,"%"PRId64"", s->n_read1+s->n_read2);
        RUN(tld_append(o,"<h2 style=\"color: red;border-bottom: none;\"> Warning: this is partial report based on the first "));
        RUN(tld_append(o,buf));

        RUN(tld_append(o," reads only.</h2>"));

        return OK;
ERROR:
        return FAIL;
}

int samstat_add_large_number(tld_strbuf *o, uint64_t n)
{
        char buf[256];
        if(n < 1000){
                snprintf(buf, 256,"%"PRId64"", n);
                RUN(tld_append(o, buf));
        }else if(n >= 1000 && n < 1000000){
                uint64_t tmp = (uint64_t) round((double)  n / 1000.0);
                snprintf(buf, 256,"%"PRId64"", tmp);
                RUN(tld_append(o, buf));
                RUN(tld_append(o, "K"));
        }else if(n >= 1000000 && n < 1000000000ull){
                uint64_t tmp = (uint64_t) round((double)  n / 1000000.0);
                snprintf(buf, 256,"%"PRId64"", tmp);
                RUN(tld_append(o, buf));
                RUN(tld_append(o, "M"));
        }else if(n >= 1000000000ull){
                uint64_t tmp = (uint64_t) round((double)  n / 1000000000.0);
                snprintf(buf, 256,"%"PRId64"", tmp);
                RUN(tld_append(o, buf));
                RUN(tld_append(o, "G"));
        }

        return OK;
ERROR:
        return FAIL;
}


int samstat_mapping_quality_overview_section(tld_strbuf *o, struct stat_collection *s)
{
        char buf[256];
        lst_node* n = NULL;

        uint64_t total = 0;
        for(int i = 0; i < 2;i++){
                for(int j = 0;j < s->mapq_map->n_bin;j++){
                        total += s->basic_nums[i][j];
                }
        }


        if(!s->is_aligned){
                RUN(tld_append(o, "<h2>File statistics</h2>\n"));
                samstat_add_large_number(o, total);
                RUN(tld_append(o, " - reads"));
                if(total >= 1000){
                        RUN(tld_append(o, " (exact: "));
                        snprintf(buf, 256,"%"PRId64"", total);
                        RUN(tld_append(o, buf));
                        RUN(tld_append(o, ")"));
                }
                RUN(tld_append(o, "<br>"));


                return OK;
        }

        /* TABLE  */
        RUN(tld_append(o, "<h2>Mapping statistics</h2>\n"));

        if(s->n_paired){
                RUN(tld_append(o, "<p>"));

                samstat_add_large_number(o, total);
                RUN(tld_append(o, " - reads "));
                if(total  >= 1000){
                        RUN(tld_append(o, " (exact: "));
                        snprintf(buf, 256,"%"PRId64"", total);
                        RUN(tld_append(o, buf));
                        RUN(tld_append(o, ")"));
                }
                RUN(tld_append(o, "<br>"));

                /* RUN(tld_append(o, " - in total<br>")); */
                if(s->n_paired){
                        /* snprintf(buf, 256,"%"PRId64"", s->n_paired); */
                        samstat_add_large_number(o, s->n_paired);
                        /* RUN(tld_append(o, buf)); */
                        RUN(tld_append(o, " - paired in sequencing ("));
                        snprintf(buf, 256,"%0.1f", (double)s->n_paired / (double) (total ) * 100.0);
                        RUN(tld_append(o, buf));
                        RUN(tld_append(o, "%)<br>"));

                        /* snprintf(buf, 256,"%"PRId64"", s->n_proper_paired); */
                        samstat_add_large_number(o, s->n_proper_paired);
                        /* RUN(tld_append(o, buf)); */
                        RUN(tld_append(o, " - properly paired ("));
                        snprintf(buf, 256,"%0.1f", (double)s->n_proper_paired / (double) s->n_paired * 100.0);
                        RUN(tld_append(o, buf));
                        RUN(tld_append(o, "%)<br>"));
                        RUN(tld_append(o, "(Stats below are based on all reads)<br>"));
                }
                RUN(tld_append(o, "</p>"));
        }

        /* LOG_MSG("%d %d -> %d  (paired: %d proper: %d )", m->n_R1_reads,m->n_R2_reads, m->n_R1_reads+m->n_R2_reads,m->n_paired, m->n_proper_paired); */


        RUN(tld_append(o, "<div id=\"MappingStatsTable\" style=\"width:100%;max-width:900px;height:240px;\"></div>\n"));
        RUN(tld_append(o,"<script>\n"));
        RUN(tld_append(o, "var mappingtableval = [\n"));
        RUN(tld_append(o, "["));
        for(int i = 0; i < s->mapq_map->n_bin;i++){
                snprintf(buf, 256, "\"%s\",",s->mapq_map->description[i]);
                RUN(tld_append(o,buf));
        }
        snprintf(buf, 256, "\"Total\",");
        RUN(tld_append(o,buf));
        o->len--;
        RUN(tld_append(o, "],"));

        total = 0;
        RUN(tld_append(o, "["));
        for(int i = 0; i < s->mapq_map->n_bin;i++){
                total += s->basic_nums[0][i];
                snprintf(buf, 256, "%"PRId64",", s->basic_nums[0][i]);
                RUN(tld_append(o,buf));
        }
        snprintf(buf, 256, "%"PRId64",", total);
        RUN(tld_append(o,buf));
        o->len--;
        RUN(tld_append(o, "],"));
        RUN(tld_append(o, "["));
        for(int i = 0; i < s->mapq_map->n_bin;i++){
                snprintf(buf, 256, "%.5f,", (double) s->basic_nums[0][i] / (double)total);
                RUN(tld_append(o,buf));
        }
        snprintf(buf, 256, "%.5f,", (double) total / (double )total);
        RUN(tld_append(o,buf));
        o->len--;
        RUN(tld_append(o, "]\n"));

        if(s->n_paired){
                RUN(tld_append(o, ","));
                uint64_t total = 0;
                RUN(tld_append(o, "["));
                for(int i = 0; i < s->mapq_map->n_bin;i++){
                        total += s->basic_nums[1][i];
                        snprintf(buf, 256, "%"PRId64",",  s->basic_nums[1][i]);
                        RUN(tld_append(o,buf));
                }
                snprintf(buf, 256, "%"PRId64",", total);
                RUN(tld_append(o,buf));
                o->len--;
                RUN(tld_append(o, "],"));
                RUN(tld_append(o, "["));
                for(int i = 0; i < s->mapq_map->n_bin;i++){
                        snprintf(buf, 256, "%.5f,", (double) s->basic_nums[1][i] / (double)total);
                        RUN(tld_append(o,buf));
                }
                snprintf(buf, 256, "%.5f,", (double) total / (double )total);
                RUN(tld_append(o,buf));
                o->len--;
                RUN(tld_append(o, "]\n"));
        }

        RUN(tld_append(o, "]\n"));

        RUN(tld_append(o, "var mappingtable_data = [{\n"));
        RUN(tld_append(o, "type: 'table',\n"));
        RUN(tld_append(o, "header: {\n"));
        if(s->n_paired){
                RUN(tld_append(o, "values: [[\"<b>Mapping Quality</b>\"], [\"<b>Read 1</b>\"], [\"<b>Percentage</b>\"], [\"<b>Read 2</b>\"], [\"<b>Percentage</b>\"]],\n"));
        }else{
                RUN(tld_append(o, "values: [[\"<b>Mapping Quality</b>\"], [\"<b>Sequences</b>\"], [\"<b>Percentage</b>\"]],\n"));
        }
        RUN(tld_append(o, "align: \"center\",\n"));
        RUN(tld_append(o, "line: {width: 1, color: 'black'},\n"));
        RUN(tld_append(o, "fill: {color: \"grey\"},\n"));
        RUN(tld_append(o, "font: {size: 14, color: \"white\"}\n"));
        RUN(tld_append(o, "},\n"));
        RUN(tld_append(o, "cells: {\n"));
        RUN(tld_append(o, "height: 40,\n"));
        RUN(tld_append(o, "values: mappingtableval,\n"));
        if(s->n_paired){
                RUN(tld_append(o, "format : [\"\", \",\", \".1%\", \",\", \".1%\"],\n"));

        }else{
                RUN(tld_append(o, "format : [\"\",\",\",\".1%\"],\n"));
        }

        RUN(tld_append(o, "align: [\"center\",\"right\"],\n"));
        RUN(tld_append(o, "line: {color: \"black\", width: 1},\n"));
        RUN(tld_append(o, "font: {size: 14, color: [\"black\"]}\n"));
        RUN(tld_append(o, "}\n"));
        RUN(tld_append(o, "}]\n"));

        RUN(tld_append(o, "var mapping_table_layout = {\n"));
        RUN(tld_append(o, "margin: {\"t\": 0, \"b\": 0, \"l\": 0, \"r\": 0},\n"));
        RUN(tld_append(o, "showlegend: false\n"));
        RUN(tld_append(o, "}\n"));


        RUN(tld_append(o, "Plotly.newPlot('MappingStatsTable', mappingtable_data,mapping_table_layout\n"));


        /* config section - here the svg buttom */



        lst_append(&n,"format : 'svg'");
        snprintf(buf, 256, "filename : 'MappingOverview'");

        lst_append(&n,buf);
        lst_append(&n,"height : 500");
        lst_append(&n,"width : 700");
        lst_append(&n, "scale : 1");


        lst_concatenate(n, ",\n");
        tld_prepend(n->data, "toImageButtonOptions : {\n");
        tld_append(n->data, "}");
        /* encapsulate whole config */
        tld_prepend(n->data, "{");
        tld_append(n->data, "}");


        RUN(tld_append(o, ",\n"));



        RUN(tld_append(o, TLD_STR(n->data)));

        lst_free(n);
        RUN(tld_append(o, ")\n"));

        RUN(tld_append(o, "</script>\n"));

        return OK;
ERROR:
        return FAIL;
}

int samstat_footer(tld_strbuf *out_buffer)
{
        RUN(tld_append(out_buffer, "<footer>\n"));

        RUN(tld_append(out_buffer, "<div>\n"));
        RUN(tld_append(out_buffer, "<section id=\"about\">\n"));
        RUN(tld_append(out_buffer, "<h3>About</h3>\n"));
        RUN(tld_append(out_buffer, "<p>SAMStat is a basic quality control program for NGS data. It displays information about mapped and unmapped reads to diagnose mapping and other problems.</p> <p> The plots on this page were draw using <a href=\"https://plotly.com/graphing-libraries/\">Plotly</a>.</p>\n"));
        RUN(tld_append(out_buffer, "</section>\n"));
        RUN(tld_append(out_buffer, "<section id=\"contact\">\n"));

        RUN(tld_append(out_buffer, "<h3>Contact</h3>\n"));

        RUN(tld_append(out_buffer, "<p>SAMStat was developed by Timo Lassmann.</p>\n"));
        RUN(tld_append(out_buffer, "<h3>Please cite:</h3>\n"));

        RUN(tld_append(out_buffer, "<p>Lassmann et al. (2011) \"SAMStat: monitoring biases in next generation sequencing data.\" Bioinformatics doi:10.1093/bioinformatics/btq614 [<a href =\"http://www.ncbi.nlm.nih.gov/pubmed/21088025/\">PMID: 21088025</a>] </p>\n"));
        RUN(tld_append(out_buffer, "</section>\n"));
        RUN(tld_append(out_buffer, "<section id=\"blogroll\">\n"));
        RUN(tld_append(out_buffer, "<h3>Links</h3>\n"));
        RUN(tld_append(out_buffer, "<ul>\n"));
        /* RUN(tld_append(out_buffer, "<li><a href=\"http://telethonkids.org.au/\">Telethon Kids Institute</a></li>\n")); */
        RUN(tld_append(out_buffer, "<li><a href=\"https://samtools.github.io/hts-specs/\">SAM/BAM specifications</a></li>\n"));
        RUN(tld_append(out_buffer, "<li><a href=\"https://github.com/TimoLassmann/SAMStat/\">SAMStat</a></li>\n"));
        /* RUN(tld_append(out_buffer, "<li><a href=\"http://code.google.com/p/bedtools/\">BEDtools</a></li>\n")); */
        RUN(tld_append(out_buffer, "<li><a href=\"https://github.com/TimoLassmann/kalign/\">Kalign</a></li>\n"));
        /* RUN(tld_append(out_buffer, "<li><a href=\"http://ngsview.sourceforge.net/\">NGSview</a></li>\n")); */


        RUN(tld_append(out_buffer, "</ul>\n"));
        RUN(tld_append(out_buffer, "</section>\n"));

        RUN(tld_append(out_buffer, "</div>\n"));
        RUN(tld_append(out_buffer, "</footer>\n"));
        return OK;
ERROR:
        return FAIL;
}
