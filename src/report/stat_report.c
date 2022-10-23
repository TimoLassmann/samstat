#include "tld.h"

#include "plot.h"
#include "../param/param.h"
#include "../collect/collect.h"
#define  STAT_REPORT_IMPORT
#include "stat_report.h"

static int samstat_footer(tld_strbuf *out_buffer);
static int write_buffer(tld_strbuf *out_buffer, struct samstat_param *p, int id);
int stat_report(struct stat_collection* s, struct samstat_param *p, int id)
{

        char* n = NULL;
        tld_strbuf* out = NULL;
        RUN(tlfilename(p->infile[id], &n));

        RUN(stat_collection_finalise(s));

        RUN(tld_strbuf_alloc(&out, 1024));
        html_header(out,n);
        plot_add(out,s->len_dist_R1);
        /* if(s->n_paired){ */
        /*          plot_add(out,s->len_dist_R2); */
        /* } */
        plot_add(out,s->base_comp_R1);
        if(s->n_paired){
                plot_add(out,s->base_comp_R2);
        }

        plot_add(out,s->qual_comp_R1);
        if(s->n_paired){
                plot_add(out,s->qual_comp_R2);
        }

        plot_add(out,s->mis_R1);
        plot_add(out,s->ins_R1);
        plot_add(out,s->del_R1);

        if(s->n_paired){
                plot_add(out,s->mis_R2);

                plot_add(out,s->ins_R2);

                plot_add(out,s->del_R2);
        }

        RUN(samstat_footer(out));
        RUN(html_end(out));

        RUN(write_buffer(out, p,id));

        tld_strbuf_free(out);
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
        RUN(tld_append(out_buffer, "<li><a href=\"https://github.com/TimoLassmann/SAMStat/\">samstat</a></li>\n"));
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
