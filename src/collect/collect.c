#include "tld.h"

#include "plot.h"
#include "sam.h"
#include "../htsinterface/htsglue.h"

/* #include "../metrics/metrics.h" */


/* #include "../sambamparse/sam_bam_parse.h" */
#define COLLECT_IMPORT
#include "collect.h"

struct mapqual_bins {
        uint8_t* map;
        char** description;
        int len;
        int n_bin;
};

#define MAPQUALBIN_gt30 0
#define MAPQUALBIN_lt30 1
#define MAPQUALBIN_lt10 2
#define MAPQUALBIN_UNMAP 3

static int get_mapqual_bins(struct mapqual_bins **map);
static void free_mapqual_bins(struct mapqual_bins *m);

/* static int collect_err(struct stat_collection* s, struct tl_seq *itm, int q_idx); */

int collect_stats(struct tl_seq_buffer *sb, struct stat_collection *s)
{
        uint64_t** t= NULL;
        int32_t q_idx;
        int32_t idx;
        /* LEt's see if we have to re-allocate stuff  */
#include "convert_tables.h"

        RUN(plot_data_resize_len(s->base_comp_R1, sb->max_len));
        RUN(plot_data_resize_len(s->base_comp_R2, sb->max_len));

        RUN(plot_data_resize_len(s->qual_comp_R1, sb->max_len));
        RUN(plot_data_resize_len(s->qual_comp_R2, sb->max_len));

        RUN(plot_data_resize_len(s->mis_R1, sb->max_len));
        RUN(plot_data_resize_len(s->mis_R2, sb->max_len));

        RUN(plot_data_resize_len(s->ins_R1, sb->max_len));
        RUN(plot_data_resize_len(s->ins_R2, sb->max_len));

        RUN(plot_data_resize_len(s->del_R1, sb->max_len));
        RUN(plot_data_resize_len(s->del_R2, sb->max_len));

        RUN(plot_data_resize_len(s->len_dist_R1, sb->max_len));
        RUN(plot_data_resize_len(s->len_dist_R2, sb->max_len));

        for(int i = 0; i < sb->num_seq;i++){
                uint8_t* seq = NULL;
                struct tl_seq *itm = sb->sequences[i];
                q_idx = MAPQUALBIN_UNMAP;
                int read = 1;
                int clip_start = 0;
                int clip_end = 0;
                if(itm->data){
                        struct aln_data* a = NULL;
                        a = itm->data;

                        if(a->flag & BAM_FPAIRED){
                                s->n_paired++;
                                if(a->flag & BAM_FPROPER_PAIR){
                                        s->n_proper_paired++;
                                }
                        }
                        if(a->reverse){
                                LOG_MSG("Oh no I should reverse!");
                                rev_comp_tl_seq(itm);
                                if(a->aln_len){
                                        reverse_comp(a->genome, a->aln_len);
                                        reverse_comp(a->read, a->aln_len);
                                }
                                a->reverse = 0;
                        }

                        q_idx = s->mapq_map->map[ a->map_q];
                        if(a->flag  & BAM_FREAD1){
                                read = 1;
                                s->n_read1++;
                        }else if(a->flag  & BAM_FREAD2){
                                read = 2;
                                s->n_read2++;
                        }
                        if(a->flag & BAM_FUNMAP){
                                q_idx = MAPQUALBIN_UNMAP;
                        }
                        clip_start = a->n_clip5;
                        clip_end = a->n_clip3;
                        if(a->aln_len){
                                uint8_t* g = NULL;
                                uint8_t* r = NULL;
                                int aln_len;
                                uint64_t** ins = NULL;
                                uint64_t** del = NULL;
                                uint64_t** mis  = NULL;

                                ins = s->ins_R1->data + s->ins_R1->group_size * q_idx;
                                del = s->del_R1->data + s->del_R1->group_size * q_idx;
                                mis = s->mis_R1->data + s->mis_R1->group_size * q_idx;
                                if(read == 2){
                                        ins = s->ins_R2->data + s->ins_R2->group_size * q_idx;
                                        del = s->del_R2->data + s->del_R2->group_size * q_idx;
                                        mis = s->mis_R2->data + s->mis_R2->group_size * q_idx;
                                }
                                aln_len = a->aln_len;
                                r = a->read;
                                g = a->genome;
                                /* fprintf(stderr,"%s\n%s\n\n",r,g); */
                                /* int m_len = MACRO_MIN((int)m->report_max_len, aln_len); */
                                int rp = 0;
                                for(int j = 0;j < aln_len;j++){
                                        /* three possibilities */

                                        /* if(r[i] == '-' && g[i] == '-'){ */

                                        /* }else  */
                                        if(r[j] != '-' && g[j] == '-'){
                                                if(j){
                                                        if(g[j-1] != '-'){
                                                                ins[0][rp]++;
                                                                /* s->ins_R1 */
                                                                /* c->ins[rp]++; */
                                                                /* c->n_ins++; */
                                                        }
                                                }else{
                                                        ins[0][rp]++;
                                                        /* c->ins[rp]++; */
                                                        /* c->n_ins++; */
                                                }
                                                rp++;
                                        }else if(r[j] == '-' && g[j] != '-'){
                                                if(j){
                                                        if(r[j-1] != '-'){
                                                                del[0][rp]++;
                                                                /* c->del[rp]++; */
                                                                /* c->n_del++; */
                                                        }
                                                }else{
                                                        del[0][rp]++;
                                                        /* c->del[rp]++; */
                                                        /* c->n_del++; */
                                                }
                                        }else if(r[j] != '-' && g[j] != '-'){
                                                if(r[j] != g[j]){
                                                        /* mis[nuc_to_internal[r[i]]][rp]++; */

                                                        switch ( r[j]) {
                                                        case 'A':
                                                        case 'a':
                                                                mis[0][rp]++;
                                                                /* c->mis[0][rp]++; */
                                                                break;
                                                        case 'C':
                                                        case 'c':
                                                                mis[1][rp]++;
                                                                /* c->mis[1][rp]++; */

                                                                break;
                                                        case 'G':
                                                        case 'g':
                                                                mis[2][rp]++;
                                                                /* c->mis[2][rp]++; */
                                                                break;
                                                        case 'T':
                                                        case 't':
                                                                mis[3][rp]++;
                                                                /* c->mis[3][rp]++; */
                                                                break;
                                                        default:

                                                                mis[4][rp]++;
                                                                /* c->mis[4][rp]++; */
                                                                break;
                                                        }
                                                        /* c->n_mis++; */
                                                }
                                                rp++;
                                        }
                                }
                        }

                }
                seq = itm->seq->str;
                /* start collecting...  */
                if(read == 1){
                        t = s->base_comp_R1->data  + s->base_comp_R1->group_size * q_idx;
                }else{
                        t = s->base_comp_R2->data  + s->base_comp_R2->group_size * q_idx;
                }
                idx = 0;
                int len = itm->len - clip_end;
                for(int i = clip_start;i < len;i++){
                        t[nuc_to_internal[seq[i]]][idx]++;
                        idx++;
                }
                /* Quality  */
                if(itm->qual->len){
                        if(itm->qual->str[0] != 0xFF){
                                if(read == 1){
                                        t = s->qual_comp_R1->data  + s->qual_comp_R1->group_size * q_idx;
                                }else{
                                        t = s->qual_comp_R2->data  + s->qual_comp_R2->group_size * q_idx;
                                }
                                seq = itm->qual->str;
                                idx = 0;
                                for(int i = clip_start;i < len;i++){
                                        int q = (int)seq[i] - sb->base_quality_offset;
                                        if(q > 41){
                                                q = 41;
                                        }
                                        /* LOG_MSG("%d %d", q, idx); */
                                        t[q][idx]++;
                                        /* data[idx][q]++; */
                                        idx++;
                                }
                        }
                }

        }
        return OK;
ERROR:
        return FAIL;
}




int stat_collection_finalise(struct stat_collection *s)
{
        ASSERT(s != NULL,"No stat collection");

        if(s->n_read2){
                RUN(tld_append(s->base_comp_R1->title," Read 1"));;
                RUN(tld_append(s->base_comp_R2->title," Read 2"));

                RUN(tld_append(s->qual_comp_R1->title," Read 1"));;
                RUN(tld_append(s->qual_comp_R2->title," Read 2"));

                RUN(tld_append(s->mis_R1->title," Read 1"));;
                RUN(tld_append(s->mis_R2->title," Read 2"));

                RUN(tld_append(s->ins_R1->title," Read 1"));;
                RUN(tld_append(s->ins_R2->title," Read 2"));

                RUN(tld_append(s->del_R1->title," Read 1"));;
                RUN(tld_append(s->del_R2->title," Read 2"));

                RUN(tld_append(s->len_dist_R1->title," Read 1"));;
                RUN(tld_append(s->len_dist_R2->title," Read 2"));
        }

        return OK;
ERROR:
        return FAIL;
}



int stat_collection_alloc(struct stat_collection **stats)
{
        struct stat_collection* s = NULL;

        MMALLOC(s,sizeof(struct stat_collection)) ;

        s->base_comp_R1 = NULL;
        s->base_comp_R2 = NULL;

        s->qual_comp_R1 = NULL;
        s->qual_comp_R2 = NULL;

        s->mis_R1 = NULL;
        s->mis_R2 = NULL;

        s->ins_R1 = NULL;
        s->ins_R2 = NULL;

        s->del_R1 = NULL;
        s->del_R2 = NULL;

        s->len_dist_R1 = NULL;
        s->len_dist_R2 = NULL;

        s->mapq_map = NULL;

        s->n_paired = 0;
        s->n_proper_paired = 0;

        s->n_read1 = 0;
        s->n_read2 = 0;


        RUN(get_mapqual_bins(&s->mapq_map));

        /* Alloc 250 len and 5 slots (ACTGN) for every mapq bin  */
        RUN(plot_data_alloc(&s->base_comp_R1, 250,  5* s->mapq_map->n_bin));
        RUN(plot_data_alloc(&s->base_comp_R2, 250,  5* s->mapq_map->n_bin));

        /* Alloc 250 len and 42 slots (base qualities ) for every mapq bin  */
        RUN(plot_data_alloc(&s->qual_comp_R1, 250,  42* s->mapq_map->n_bin));
        RUN(plot_data_alloc(&s->qual_comp_R2, 250,  42* s->mapq_map->n_bin));

        /* Alloc 250 len and 5 slots (ACTGN) for every mapq bin  */
        RUN(plot_data_alloc(&s->mis_R1, 250,  5* s->mapq_map->n_bin));
        RUN(plot_data_alloc(&s->mis_R2, 250,  5* s->mapq_map->n_bin));

        /* Alloc 250 len and 1 slot for every mapq bin */
        RUN(plot_data_alloc(&s->ins_R1, 250,  s->mapq_map->n_bin));
        RUN(plot_data_alloc(&s->ins_R2, 250,  s->mapq_map->n_bin));

        /* Alloc 250 len and 1 slot for every mapq bin */
        RUN(plot_data_alloc(&s->del_R1, 250,  s->mapq_map->n_bin));
        RUN(plot_data_alloc(&s->del_R2, 250,  s->mapq_map->n_bin));

        /* Alloc 250 len and 1 slot for every mapq bin */
        RUN(plot_data_alloc(&s->len_dist_R1, 250,  s->mapq_map->n_bin));
        RUN(plot_data_alloc(&s->len_dist_R2, 250,  s->mapq_map->n_bin));

        /* Configure plots  */

        /* base composition */

        char** base_label = NULL;

        galloc(&base_label, 5,2);
        snprintf(base_label[0], 2, "%s","A");
        snprintf(base_label[1], 2, "%s","C");
        snprintf(base_label[2], 2, "%s","G");
        snprintf(base_label[3], 2, "%s","T");
        snprintf(base_label[4], 2, "%s","N");


        /* s->mapq_map->description */
        RUN(plot_data_config(s->base_comp_R1,
                             PLOT_TYPE_BAR,
                             PLOT_MOD_NORMAL,
                             5,
                             PLOT_VIZ_FIRSTGROUP,
                             "bc_R1",
                             "Base Composition",
                             "Length",
                             "Counts",
                             "BaseComposition",
                             base_label,
                             s->mapq_map->description
                    ));

        RUN(plot_data_config(s->base_comp_R2,
                             PLOT_TYPE_BAR,
                             PLOT_MOD_NORMAL,
                             5,
                             PLOT_VIZ_FIRSTGROUP,
                             "bc_R2",
                             "Base Composition",
                             "Length",
                             "Counts",
                             "BaseCompositionR2",
                             base_label,
                             s->mapq_map->description
                    ));


        /* qual composition */

        RUN(plot_data_config(s->qual_comp_R1,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_ERROR_BAR,
                             42,
                             PLOT_VIZ_ALL,
                             "bq_R1",
                             "Base quality distribution",
                             "Length",
                             "Base Quality",
                             "QualComposition",
                             NULL,
                             s->mapq_map->description
                    ));

        RUN(plot_data_config(s->qual_comp_R2 ,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_ERROR_BAR,
                             42,
                             PLOT_VIZ_ALL,
                             "bq_R2",
                             "Base quality distribution",
                             "Length",
                             "Base Quality",
                             "QualCompositionR2",
                             NULL,
                             s->mapq_map->description
                    ));


        /* mismatch composition */
        /* LOG_MSG("MIS comp"); */
        RUN(plot_data_config(s->mis_R1,
                             PLOT_TYPE_BAR,
                             PLOT_MOD_NORMAL,
                             5,
                             PLOT_VIZ_FIRSTGROUP,
                             "mis_R1",
                             "Distribution of Mismatches",
                             "Length",
                             "Counts",
                             "MismatchDistribution",
                             base_label,
                             s->mapq_map->description
                    ));

        RUN(plot_data_config(s->mis_R2,
                             PLOT_TYPE_BAR,
                             PLOT_MOD_NORMAL,
                             5,
                             PLOT_VIZ_FIRSTGROUP,
                             "mis_R2",
                             "Distribution of Mismatches",
                             "Length",
                             "Counts",
                             "MismatchDistributionR2",
                             base_label,
                             s->mapq_map->description
                    ));

        /* ins composition */
        /* LOG_MSG("Ins comp"); */
        RUN(plot_data_config(s->ins_R1,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_NORMAL,
                             1,//s->mapq_map->n_bin,
                             PLOT_VIZ_ALL,
                             "ins_R1",
                             "Distribution of Insertions",
                             "Length",
                             "Counts",
                             "InsertionDistribution",
                             s->mapq_map->description,
                             NULL
                    ));

        RUN(plot_data_config(s->ins_R2,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_NORMAL,
                             1,//s->mapq_map->n_bin,
                             PLOT_VIZ_ALL,
                             "ins_R2",
                             "Distribution of Insertions",
                             "Length",
                             "Counts",
                             "InsertionDistribution",
                             s->mapq_map->description,
                             NULL
                    ));

        /* del composition */
        /* LOG_MSG("Del comp"); */
        RUN(plot_data_config(s->del_R1,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_NORMAL,
                             1,//s->mapq_map->n_bin,
                             PLOT_VIZ_ALL,
                             "del_R1",
                             "Distribution of Deletions",
                             "Length",
                             "Counts",
                             "DeletionDistribution",
                             s->mapq_map->description,
                             NULL
                    ));

        RUN(plot_data_config(s->del_R2,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_NORMAL,
                             1,//s->mapq_map->n_bin,
                             PLOT_VIZ_ALL,
                             "del_R2",
                             "Distribution of Deletions",
                             "Length",
                             "Counts",
                             "DeletionDistribution",
                             s->mapq_map->description,
                             NULL
                    ));

        gfree(base_label);
        *stats = s;

        return OK;
ERROR:
        return FAIL;
}

void stat_collection_free(struct stat_collection *s)
{
        if(s){
                if(s->base_comp_R1){
                        plot_data_free(s->base_comp_R1);
                }
                if(s->base_comp_R2){
                        plot_data_free(s->base_comp_R2);
                }

                if(s->qual_comp_R1){
                        plot_data_free(s->qual_comp_R1);
                }
                if(s->qual_comp_R2){
                        plot_data_free(s->qual_comp_R2);
                }

                if(s->mis_R1){
                        plot_data_free(s->mis_R1);
                }
                if(s->mis_R2){
                        plot_data_free(s->mis_R2);
                }

                if(s->ins_R1){
                        plot_data_free(s->ins_R1);
                }
                if(s->ins_R2){
                        plot_data_free(s->ins_R2);
                }

                if(s->del_R1){
                        plot_data_free(s->del_R1);
                }
                if(s->del_R2){
                        plot_data_free(s->del_R2);
                }

                if(s->len_dist_R1){
                        plot_data_free(s->len_dist_R1);
                }
                if(s->len_dist_R2){
                        plot_data_free(s->len_dist_R2);
                }

                if(s->mapq_map){
                        free_mapqual_bins(s->mapq_map);
                }
                MFREE(s);
        }
}


int get_mapqual_bins(struct mapqual_bins **map)
{
        struct mapqual_bins* m = NULL;

        MMALLOC(m, sizeof(struct mapqual_bins));

        m->map = NULL;
        m->description = NULL;
        m->len = 256;
        m->n_bin = 4;
        galloc(&m->map,m->len);
        galloc(&m->description, m->n_bin, 256);

        for(int i = 0; i < m->len;i++){
                /* if(i == 0){ */
                /*         m->map[i] = MAPQUALBIN_UNMAP; */
                /* }else  */if(i < 10){
                        m->map[i] = MAPQUALBIN_lt10;
                }else if(i < 30){
                        m->map[i] = MAPQUALBIN_lt30;
                }else{
                        m->map[i] = MAPQUALBIN_gt30;
                }
        }

        snprintf(m->description[MAPQUALBIN_gt30], 256, "MAPQ gt 30");
        snprintf(m->description[MAPQUALBIN_lt30], 256, "MAPQ lt 30");
        snprintf(m->description[MAPQUALBIN_lt10], 256, "MAPQ lt 10");
        snprintf(m->description[MAPQUALBIN_UNMAP], 256, "Unmapped");


        *map = m;
        return OK;
ERROR:
        free_mapqual_bins(m);
        return FAIL;
}

void free_mapqual_bins(struct mapqual_bins *m)
{
        if(m){
                gfree(m->description);
                gfree(m->map);
                MFREE(m);
        }

}
