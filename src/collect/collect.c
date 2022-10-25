#include "tld.h"

#include "plot.h"
#include "sam.h"
#include "../htsinterface/htsglue.h"
#include "../param/param.h"

#include "convert_tables.h"

/* #include "convert_tables.h" */
/* #include "../metrics/metrics.h" */


/* #include "../sambamparse/sam_bam_parse.h" */
#define COLLECT_IMPORT
#include "collect.h"


#define MAPQUALBIN_gt30 0
#define MAPQUALBIN_lt30 1
#define MAPQUALBIN_lt10 2
#define MAPQUALBIN_UNMAP 3

#define COLLECT_TYPE_SEQ 1
#define COLLECT_TYPE_SEQQUAL 2
#define COLLECT_TYPE_ALL 3

static int get_mapqual_bins(struct mapqual_bins **map);
static void free_mapqual_bins(struct mapqual_bins *m);

static int config_stats(struct tl_seq_buffer *sb, struct stat_collection *s);
static inline int get_aln_errors(struct stat_collection *s, struct aln_data *a,
                                 int read, int q_idx);


static int collect_stats_seq(struct tl_seq_buffer *sb, struct stat_collection *s);
static int collect_stats_seqqual(struct tl_seq_buffer *sb, struct stat_collection *s);
static int collect_stats_all(struct tl_seq_buffer *sb, struct stat_collection *s);

static int collect_stats_end_seqqual(struct tl_seq_buffer *sb, struct stat_collection *s);


int collect_stats(struct tl_seq_buffer *sb, struct stat_collection *s)
{
        /* uint64_t** t= NULL; */
        /* int32_t q_idx; */
        /* int32_t idx; */
        /* LEt's see if we have to re-allocate stuff  */
        if(!s->config){
                config_stats(sb, s);
        }


        RUN(plot_data_resize_len(s->base_comp_R1, sb->max_len+1));
        RUN(plot_data_resize_len(s->base_comp_R2, sb->max_len+1));

        RUN(plot_data_resize_len(s->qual_comp_R1, sb->max_len+1));
        RUN(plot_data_resize_len(s->qual_comp_R2, sb->max_len+1));

        RUN(plot_data_resize_len(s->mis_R1, sb->max_len+1));
        RUN(plot_data_resize_len(s->mis_R2, sb->max_len+1));

        RUN(plot_data_resize_len(s->ins_R1, sb->max_len+1));
        RUN(plot_data_resize_len(s->ins_R2, sb->max_len+1));

        RUN(plot_data_resize_len(s->del_R1, sb->max_len+1));
        RUN(plot_data_resize_len(s->del_R2, sb->max_len+1));

        RUN(plot_data_resize_len(s->len_dist_R1, sb->max_len+1));
        RUN(plot_data_resize_len(s->len_dist_R2, sb->max_len+1));

        if(s->collect_end){
                RUN(plot_data_resize_len(s->base_comp_R1_end, sb->max_len+1));
                RUN(plot_data_resize_len(s->base_comp_R2_end, sb->max_len+1));

                RUN(plot_data_resize_len(s->qual_comp_R1_end, sb->max_len+1));
                RUN(plot_data_resize_len(s->qual_comp_R2_end, sb->max_len+1));

        }

        
        /* additional plots will not be resized...  */
        /* collect_stats_all(sb,s); */

        switch (s->collect_type) {
        case COLLECT_TYPE_ALL:
                collect_stats_all(sb,s);
                if(s->collect_end){
                        collect_stats_end_seqqual(sb,s);
                }
                break;
        case COLLECT_TYPE_SEQQUAL:
                collect_stats_seqqual(sb,s);
                if(s->collect_end){
                        collect_stats_end_seqqual(sb,s);
                }
                break;
        case COLLECT_TYPE_SEQ:
                collect_stats_seq(sb,s);
                break;
        default:
                collect_stats_all(sb,s);
        }

        /* for(int i = 0; i < sb->num_seq;i++){ */
        /*         uint8_t* seq = NULL; */
        /*         struct tl_seq *itm = sb->sequences[i]; */
        /*         q_idx = MAPQUALBIN_UNMAP; */
        /*         int read = 1; */
        /*         int clip_start = 0; */
        /*         int clip_end = 0; */
        /*         if(itm->data){ */
        /*                 struct aln_data* a = NULL; */
        /*                 a = itm->data; */

        /*                 if(a->flag & BAM_FPAIRED){ */
        /*                         s->n_paired++; */
        /*                         if(a->flag & BAM_FPROPER_PAIR){ */
        /*                                 s->n_proper_paired++; */
        /*                         } */
        /*                 } */
        /*                 if(a->reverse){ */
        /*                         LOG_MSG("Oh no I should reverse!"); */
        /*                         rev_comp_tl_seq(itm); */
        /*                         if(a->aln_len){ */
        /*                                 reverse_comp(a->genome, a->aln_len); */
        /*                                 reverse_comp(a->read, a->aln_len); */
        /*                         } */
        /*                         a->reverse = 0; */
        /*                 } */

        /*                 q_idx = s->mapq_map->map[ a->map_q]; */
        /*                 if(a->flag  & BAM_FREAD1){ */
        /*                         read = 1; */
        /*                         s->n_read1++; */
        /*                 }else if(a->flag  & BAM_FREAD2){ */
        /*                         read = 2; */
        /*                         s->n_read2++; */
        /*                 } */
        /*                 if(a->flag & BAM_FUNMAP){ */
        /*                         q_idx = MAPQUALBIN_UNMAP; */
        /*                 } */
        /*                 clip_start = a->n_clip5; */
        /*                 clip_end = a->n_clip3; */
        /*                 if(a->aln_len){ */
        /*                         get_aln_errors(s, a, read, q_idx); */
        /*                 } */
        /*         } */
        /*         if(read == 1){ */
        /*                 t = s->len_dist_R1->data;// + q_idx; */
        /*                 s->basic_nums[0][q_idx]++; */
        /*         }else{ */
        /*                 t = s->len_dist_R2->data;// +  q_idx; */
        /*                 s->basic_nums[1][q_idx]++; */
        /*         } */
        /*         /\* LOG_MSG("Adding : %d %d %d  ", itm->len - (clip_start + clip_end),  s->len_dist_R1->len,s->len_dist_R2->len ); *\/ */
        /*         t[q_idx][itm->len - (clip_start + clip_end)]++; */

        /*         seq = itm->seq->str; */


        /*         /\* start collecting...  *\/ */
        /*         if(read == 1){ */
        /*                 t = s->base_comp_R1->data  + s->base_comp_R1->group_size * q_idx; */
        /*         }else{ */
        /*                 t = s->base_comp_R2->data  + s->base_comp_R2->group_size * q_idx; */
        /*         } */
        /*         idx = 0; */
        /*         int len = itm->len - clip_end; */
        /*         for(int i = clip_start;i < len;i++){ */
        /*                 t[nuc_to_internal[seq[i]]][idx]++; */
        /*                 idx++; */
        /*         } */
        /*         /\* Quality  *\/ */
        /*         if(itm->qual->len){ */
        /*                 if(itm->qual->str[0] != 0xFF){ */
        /*                         if(read == 1){ */
        /*                                 t = s->qual_comp_R1->data  + s->qual_comp_R1->group_size * q_idx; */
        /*                         }else{ */
        /*                                 t = s->qual_comp_R2->data  + s->qual_comp_R2->group_size * q_idx; */
        /*                         } */
        /*                         seq = itm->qual->str; */
        /*                         idx = 0; */
        /*                         for(int i = clip_start;i < len;i++){ */
        /*                                 int q = (int)seq[i] - sb->base_quality_offset; */
        /*                                 if(q > 41){ */
        /*                                         q = 41; */
        /*                                 } */
        /*                                 /\* LOG_MSG("%d %d", q, idx); *\/ */
        /*                                 t[q][idx]++; */
        /*                                 /\* data[idx][q]++; *\/ */
        /*                                 idx++; */
        /*                         } */
        /*                 } */
        /*         } */
        /* } */
        return OK;
ERROR:
        return FAIL;
}

int config_stats(struct tl_seq_buffer *sb, struct stat_collection *s)
{
        uint8_t has_aln_data = 0;
        uint8_t has_qual_data = 0;

        for(int i = 0; i < sb->num_seq;i++){
                struct tl_seq *itm = sb->sequences[i];
                if(itm->data){
                        has_aln_data = 1;
                }
                if(itm->qual->len){
                        if(itm->qual->str[0] != 0xFF){
                                has_qual_data = 1;
                        }
                }

        }
        if(has_aln_data && has_qual_data ){
                s->collect_type = COLLECT_TYPE_ALL;
        }else if(!has_aln_data && has_qual_data){
                s->collect_type = COLLECT_TYPE_SEQQUAL;
        }else if(!has_aln_data && !has_qual_data){
                s->collect_type = COLLECT_TYPE_SEQ;
        }else{
                s->collect_type = COLLECT_TYPE_ALL;
        }

        s->config = 1;
        return OK;
}

int collect_stats_seq(struct tl_seq_buffer *sb, struct stat_collection *s)
{
        uint64_t** t= NULL;
        int32_t q_idx;
        q_idx = MAPQUALBIN_UNMAP;

        for(int i = 0; i < sb->num_seq;i++){
                uint8_t* seq = NULL;

                struct tl_seq *itm = sb->sequences[i];


                t = s->len_dist_R1->data;// + q_idx;
                s->basic_nums[0][q_idx]++;

                /* LOG_MSG("Adding : %d %d %d  ", itm->len - (clip_start + clip_end),  s->len_dist_R1->len,s->len_dist_R2->len ); */
                t[q_idx][itm->len]++;

                seq = itm->seq->str;


                /* start collecting...  */
                t = s->base_comp_R1->data  + s->base_comp_R1->group_size * q_idx;
                int len = itm->len;
                for(int i = 0;i < len;i++){
                        t[nuc_to_internal[seq[i]]][i]++;
                }
        }
        return OK;
}


int collect_stats_seqqual(struct tl_seq_buffer *sb, struct stat_collection *s)
{
        uint64_t** t= NULL;
        uint64_t** t_q= NULL;
        int32_t q_idx;

        q_idx = MAPQUALBIN_UNMAP;
        for(int i = 0; i < sb->num_seq;i++){
                uint8_t* seq = NULL;
                uint8_t* qual = NULL;
                struct tl_seq *itm = sb->sequences[i];

                t = s->len_dist_R1->data;// + q_idx;
                s->basic_nums[0][q_idx]++;

                /* LOG_MSG("Adding : %d %d %d  ", itm->len - (clip_start + clip_end),  s->len_dist_R1->len,s->len_dist_R2->len ); */
                t[q_idx][itm->len]++;

                seq = itm->seq->str;
                qual = itm->qual->str;

                /* start collecting...  */
                t = s->base_comp_R1->data  + s->base_comp_R1->group_size * q_idx;
                t_q = s->qual_comp_R1->data  + s->qual_comp_R1->group_size * q_idx;
                int len = itm->len;
                for(int i = 0;i < len;i++){
                        t[nuc_to_internal[seq[i]]][i]++;
                        int q = (int)qual[i] - sb->base_quality_offset;
                        if(q > 41){
                                q = 41;
                        }
                        t_q[q][i]++;
                }
        }
        return OK;
}


int collect_stats_all(struct tl_seq_buffer *sb, struct stat_collection *s)
{
        uint64_t** t= NULL;
        int32_t q_idx;
        int32_t idx;
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
                                get_aln_errors(s, a, read, q_idx);
                        }
                }
                if(read == 1){
                        t = s->len_dist_R1->data;// + q_idx;
                        s->basic_nums[0][q_idx]++;
                }else{
                        t = s->len_dist_R2->data;// +  q_idx;
                        s->basic_nums[1][q_idx]++;
                }
                /* LOG_MSG("Adding : %d %d %d  ", itm->len - (clip_start + clip_end),  s->len_dist_R1->len,s->len_dist_R2->len ); */
                t[q_idx][itm->len - (clip_start + clip_end)]++;

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
}

int collect_stats_end_seqqual(struct tl_seq_buffer *sb, struct stat_collection *s)
{
        uint64_t** t= NULL;
        uint64_t** t_q= NULL;
        int32_t q_idx;

        q_idx = MAPQUALBIN_UNMAP;
        for(int i = 0; i < sb->num_seq;i++){
                uint8_t* seq = NULL;
                uint8_t* qual = NULL;
                struct tl_seq *itm = sb->sequences[i];

                q_idx = MAPQUALBIN_UNMAP;
                int read = 1;
                if(itm->data){
                        struct aln_data* a = NULL;
                        a = itm->data;
                        q_idx = s->mapq_map->map[ a->map_q];
                        if(a->flag  & BAM_FREAD1){
                                read = 1;
                        }else if(a->flag  & BAM_FREAD2){
                                read = 2;
                        }
                        if(a->flag & BAM_FUNMAP){
                                q_idx = MAPQUALBIN_UNMAP;
                        }
                }
                seq = itm->seq->str;
                qual = itm->qual->str;

                /* start collecting...  */
                if(read == 1){
                        t = s->base_comp_R1_end->data  + s->base_comp_R1->group_size * q_idx;
                        t_q = s->qual_comp_R1_end->data  + s->qual_comp_R1->group_size * q_idx;
                }else{
                        t = s->base_comp_R2_end->data  + s->base_comp_R2->group_size * q_idx;
                        t_q = s->qual_comp_R2_end->data  + s->qual_comp_R2->group_size * q_idx;
                }

                ASSERT(s->base_comp_R1_end->len ==s->qual_comp_R1_end->len,"Lengths differ seq / qual");
                int plot_len = s->base_comp_R1_end->len;
                int len = itm->len-1;
                int idx = 0;

                while(len >= 0 && idx < plot_len){
                        /* LOG_MSG("idx = %d", idx); */
                        t[nuc_to_internal[seq[len]]][idx]++;
                        int q = (int)qual[len] - sb->base_quality_offset;
                        if(q > 41){
                                q = 41;
                        }
                        t_q[q][idx]++;
                        idx++;
                        len--;
                }
        }
        return OK;
ERROR:
        return FAIL;
}


inline int get_aln_errors(struct stat_collection* s,struct aln_data* a, int read,int q_idx)
{
        uint8_t* g = NULL;
        uint8_t* r = NULL;
        int aln_len;
        uint64_t** ins = NULL;
        uint64_t** del = NULL;
        uint64_t** mis  = NULL;


        ins = s->ins_R1->data;// + s->ins_R1->group_size * q_idx;
        del = s->del_R1->data;// + s->del_R1->group_size * q_idx;
        mis = s->mis_R1->data + s->mis_R1->group_size * q_idx;
        if(read == 2){
                ins = s->ins_R2->data;// + s->ins_R2->group_size * q_idx;
                del = s->del_R2->data;// + s->del_R2->group_size * q_idx;
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
                                        ins[q_idx][rp]++;
                                }
                        }else{
                                ins[q_idx][rp]++;
                        }
                        rp++;
                }else if(r[j] == '-' && g[j] != '-'){
                        if(j){
                                if(r[j-1] != '-'){
                                        del[q_idx][rp]++;
                                }
                        }else{
                                del[q_idx][rp]++;
                        }
                }else if(r[j] != '-' && g[j] != '-'){
                        if(r[j] != g[j]){
                                mis[nuc_to_internal[r[j]]][rp]++;
                        }
                        rp++;
                }
        }
        return OK;
/* ERROR: */
/*         return FAIL; */
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

int stat_collection_config_additional_plot(struct stat_collection *s,struct samstat_param *p)
{
        if(p->collect_end){
                s->collect_end = 1;
                RUN(plot_data_alloc(&s->base_comp_R1_end, 250,  5* s->mapq_map->n_bin));
                RUN(plot_data_alloc(&s->base_comp_R2_end, 250,  5* s->mapq_map->n_bin));

                /* Alloc 250 len and 42 slots (base qualities ) for every mapq bin  */
                RUN(plot_data_alloc(&s->qual_comp_R1_end, 250,  42* s->mapq_map->n_bin));
                RUN(plot_data_alloc(&s->qual_comp_R2_end, 250,  42* s->mapq_map->n_bin));

                char** base_label = NULL;

                galloc(&base_label, 5,2);
                snprintf(base_label[0], 2, "%s","A");
                snprintf(base_label[1], 2, "%s","C");
                snprintf(base_label[2], 2, "%s","G");
                snprintf(base_label[3], 2, "%s","T");
                snprintf(base_label[4], 2, "%s","N");


                /* s->mapq_map->description */
                RUN(plot_data_config(s->base_comp_R1_end,
                                     PLOT_TYPE_BAR,
                                     PLOT_MOD_NORMAL,
                                     5,
                                     PLOT_EMPTY_SERIES, /* plot empty series - yes please  */
                                     "bc_R1_end",
                                     "Base Composition at read end",
                                     "Position relative to read end",
                                     "Counts",
                                     "BaseCompositionend",
                                     base_label,
                                     s->mapq_map->description
                            ));

                RUN(plot_data_config(s->base_comp_R2_end,
                                     PLOT_TYPE_BAR,
                                     PLOT_MOD_NORMAL,
                                     5,
                                     PLOT_EMPTY_SERIES, /* plot empty series - yes please  */
                                     "bc_R2_end",
                                     "Base Composition at read end",
                                     "Position relative to read end",
                                     "Counts",
                                     "BaseCompositionR2end",
                                     base_label,
                                     s->mapq_map->description
                            ));



                /* qual composition */

                RUN(plot_data_config(s->qual_comp_R1_end,
                                     PLOT_TYPE_LINES,
                                     PLOT_MOD_ERROR_BAR,
                                     42,
                                     0,
                                     "bq_R1_end",
                                     "Base quality distribution at read end",
                                     "Position relative to read end",
                                     "Base Quality",
                                     "QualCompositionend",
                                     NULL,
                                     s->mapq_map->description
                            ));

                RUN(plot_data_config(s->qual_comp_R2_end ,
                                     PLOT_TYPE_LINES,
                                     PLOT_MOD_ERROR_BAR,
                                     42,
                                     0,
                                     "bq_R2_end",
                                     "Base quality distribution at read end",
                                     "Position relative to read end",
                                     "Base Quality",
                                     "QualCompositionR2end",
                                     NULL,
                                     s->mapq_map->description
                            ));
                gfree(base_label);
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

        /* optional plots  */
        s->base_comp_R1_end = NULL;
        s->base_comp_R2_end = NULL;

        s->qual_comp_R1_end = NULL;
        s->qual_comp_R2_end = NULL;

        s->mapq_map = NULL;

        s->basic_nums = NULL;

        s->n_paired = 0;
        s->n_proper_paired = 0;

        s->n_read1 = 0;
        s->n_read2 = 0;

        s->is_aligned = 0;
        s->is_partial_report = 0;
        s->config = 0;
        s->collect_type = 0;
        s->collect_end = 0;
        RUN(get_mapqual_bins(&s->mapq_map));


        galloc(&s->basic_nums, 2, s->mapq_map->n_bin);
        for(int i = 0; i < 2;i++){
                for(int j = 0; j < s->mapq_map->n_bin;j++){
                        s->basic_nums[i][j] = 0;
                }
        }

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
                             PLOT_EMPTY_SERIES, /* plot empty series - yes please  */
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
                             PLOT_EMPTY_SERIES, /* plot empty series - yes please  */
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
                             0,
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
                             0,
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
                             PLOT_EMPTY_SERIES,
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
                             PLOT_EMPTY_SERIES,
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
                             0,//s->mapq_map->n_bin,
                             0,
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
                             0,//s->mapq_map->n_bin,
                             0,
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
                             0,//s->mapq_map->n_bin,
                             0,
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
                             0,//s->mapq_map->n_bin,
                             0,
                             "del_R2",
                             "Distribution of Deletions",
                             "Length",
                             "Counts",
                             "DeletionDistribution",
                             s->mapq_map->description,
                             NULL
                    ));

        /* LOG_MSG("LEN IDST "); */
        /* Length dist  */
        RUN(plot_data_config(s->len_dist_R1,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_DENSITY,
                             0,//s->mapq_map->n_bin,
                             0,
                             "len_R1",
                             "Length Distribution",
                             "Length",
                             "Counts",
                             "LengthDistribution",
                             s->mapq_map->description,
                             NULL

                    ));

        RUN(plot_data_config(s->len_dist_R2,
                             PLOT_TYPE_LINES,
                             PLOT_MOD_DENSITY,
                             0,//s->mapq_map->n_bin,
                             0,
                             "len_R2",
                             "Length Distribution",
                             "Length",
                             "Counts",
                             "LengthDistribution",
                             s->mapq_map->description,
                             NULL
                    ));

        /* fprintf(stderr,"%s -> %d\n", TLD_STR(s->len_dist_R1->title), s->len_dist_R1->group_size); */
        /* fprintf(stderr,"%s -> %d\n", TLD_STR(s->base_comp_R1->title), s->base_comp_R1->group_size); */
        /* fprintf(stderr,"%s -> %d\n", TLD_STR(s->qual_comp_R1->title), s->qual_comp_R1->group_size); */
        /* fprintf(stderr,"%s -> %d\n", TLD_STR(s->mis_R1->title), s->mis_R1->group_size); */
        /* fprintf(stderr,"%s -> %d\n", TLD_STR(s->ins_R1->title), s->ins_R1->group_size); */
        /* fprintf(stderr,"%s -> %d\n", TLD_STR(s->del_R1->title), s->del_R1->group_size); */


        /* exit(0); */
        gfree(base_label);
        *stats = s;

        return OK;
ERROR:
        return FAIL;
}

void stat_collection_free(struct stat_collection *s)
{
        if(s){
                if(s->basic_nums){
                        gfree(s->basic_nums);
                }
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
                /* Additional plots */

                if(s->base_comp_R1_end){
                        plot_data_free(s->base_comp_R1_end);
                }
                if(s->base_comp_R2_end){
                        plot_data_free(s->base_comp_R2_end);
                }

                if(s->qual_comp_R1_end){
                        plot_data_free(s->qual_comp_R1_end);
                }
                if(s->qual_comp_R2_end){
                        plot_data_free(s->qual_comp_R2_end);
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
