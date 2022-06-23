

#include "H5public.h"
#include "alloc/tld-alloc.h"
#include "core/tld-core.h"
#include "sam.h"
#include "seq/tld-seq.h"
#include "tld.h"

#include "../sambamparse/sam_bam_parse.h"
#include "../htsinterface/htsglue.h"
#include <limits.h>
#include <stdint.h>
#include <stdio.h>

#define  METRICS_IMPORT
#include "metrics.h"

/* #define REPORT_MAX_LEN 500 */

static int collect_seq_comp(struct metrics *m, struct tl_seq *s);
static int collect_qual_comp(struct metrics *m, struct tl_seq *s, const int offset);
static int collect_error_comp(struct metrics *m, struct tl_seq *s);
static int collect_len_comp(struct metrics *m, struct tl_seq *s);

static int set_read_mapqbin( struct tl_seq *s, struct mapqual_bins* map, int *read, int *idx);

static int alloc_seq_comp(struct seq_composition **seq_comp, int L,int report_max_len);
/* static int resize_seq_comp(struct seq_composition* s,int newL, int new_max_len); */
static void free_seq_comp(struct seq_composition *s);

static int alloc_qual_comp(struct qual_composition **qual_comp,int report_max_len);
/* static int resize_qual_comp(struct qual_composition* q, int new_max_len); */
static void free_qual_comp(struct qual_composition *q);

static int alloc_error_comp(struct error_composition** error_comp, int L,int report_max_len);
/* static int resize_error_comp(struct error_composition *e, int newL, int new_max_len); */
static void free_error_comp(struct error_composition *e);

static int alloc_len_comp(struct len_composition** len_comp,int cur_max_len);
static int resize_len_comp(struct len_composition* s,int new_max_len);
static void free_len_comp(struct len_composition *s);

static int get_mapqual_bins(struct mapqual_bins **map);
static void free_mapqual_bins(struct mapqual_bins *m);


int get_metrics(struct tl_seq_buffer *sb, struct metrics *m)
{
        int len_change = 0;
        int read = 0;
        for(int i = 0; i < sb->num_seq;i++){
                read = 0;
                if(sb->sequences[i]->data){
                        struct aln_data* a = NULL;
                        a = sb->sequences[i]->data;
                        if(a->flag & BAM_FPAIRED){
                                m->n_paired++;
                                if(a->flag & BAM_FPROPER_PAIR){
                                        m->n_proper_paired++;
                                }
                        }

                        if(a->flag & BAM_FREAD1){
                                m->n_R1_reads++;
                                if(m->max_len_R1 < (uint32_t)sb->sequences[i]->len){
                                        m->max_len_R1 = sb->sequences[i]->len;
                                        len_change = 1;
                                }
                                if(m->min_len_R1 > (uint32_t)sb->sequences[i]->len){
                                        m->min_len_R1 = sb->sequences[i]->len;
                                }
                        }else if(a->flag & BAM_FREAD2){
                                m->n_R2_reads++;
                                if(m->max_len_R2 < (uint32_t)sb->sequences[i]->len){
                                        m->max_len_R2 = sb->sequences[i]->len;
                                        len_change = 2;
                                }
                                if(m->min_len_R2 > (uint32_t)sb->sequences[i]->len){
                                        m->min_len_R2 = sb->sequences[i]->len;
                                }
                        }else{
                                /* LOG_MSG("Huh?"); */
                                m->n_R1_reads++;
                                if(m->max_len_R1 < (uint32_t)sb->sequences[i]->len){
                                        m->max_len_R1 = sb->sequences[i]->len;
                                        len_change = 1;
                                }
                                if(m->min_len_R1 > (uint32_t)sb->sequences[i]->len){
                                        m->min_len_R1 = sb->sequences[i]->len;
                                }
                        }
                }else{
                        m->n_R1_reads++;
                        if(m->max_len_R1 < (uint32_t)sb->sequences[i]->len){
                                m->max_len_R1 = sb->sequences[i]->len;
                                len_change = 1;
                        }
                        if(m->min_len_R1 > (uint32_t)sb->sequences[i]->len){
                                m->min_len_R1 = sb->sequences[i]->len;
                        }
                }
        }
        /* LOG_MSG("Max len: %d %d", m->max_len_R1, m->max_len_R2); */
        /* Sequence composition */
        for(int i = 0; i < m->n_mapq_bins;i++){
                if(!m->seq_comp_R1[i]){
                        RUN(alloc_seq_comp(&m->seq_comp_R1[i], sb->L,m->report_max_len));
                        RUN(alloc_qual_comp(&m->qual_comp_R1[i],m->report_max_len));
                        RUN(alloc_error_comp(&m->error_comp_R1[i], sb->L,m->report_max_len));
                        RUN(alloc_len_comp(&m->len_comp_R1[i], m->max_len_R1));
                } else if(len_change == 1){
                        RUN(resize_len_comp(m->len_comp_R1[i], m->max_len_R1));
                }
                                /*         /\* LOG_MSG("Changing length: %d L %d ", m->max_len_R1, sb->L); *\/ */
                                /*         RUN(resize_seq_comp(m->seq_comp_R1[i], sb->L, m->max_len_R1)); */
                                /*         RUN(resize_qual_comp(m->qual_comp_R1[i], m->max_len_R1)); */
                                /*         RUN(resize_error_comp(m->error_comp_R1[i], sb->L, m->max_len_R1)); */

        }
        if(m->n_paired){
                for(int i = 0; i < m->n_mapq_bins;i++){
                        if(!m->seq_comp_R2[i]){
                                RUN(alloc_seq_comp(&m->seq_comp_R2[i], sb->L,m->report_max_len));
                                RUN(alloc_qual_comp(&m->qual_comp_R2[i],m->report_max_len));
                                RUN(alloc_error_comp(&m->error_comp_R2[i], sb->L,m->report_max_len));
                                RUN(alloc_len_comp(&m->len_comp_R2[i], m->max_len_R2));
                        } else if(len_change == 2){
                                RUN(resize_len_comp(m->len_comp_R2[i], m->max_len_R2));
                        }
                }
        }

        for(int i = 0; i < sb->num_seq;i++){
                collect_seq_comp(m, sb->sequences[i]);
                if(sb->sequences[i]->qual->len){
                        collect_qual_comp(m, sb->sequences[i], sb->base_quality_offset);
                }
                collect_error_comp(m, sb->sequences[i]);
                collect_len_comp(m, sb->sequences[i]);
        }


        return OK;
ERROR:
        return FAIL;
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



/* int debug_metrics_print(struct metrics *m) */
/* { */
/*         ASSERT(m != NULL,"No metrics"); */
/*         LOG_MSG("Min len : %d", m->min_len); */
/*         LOG_MSG("Max len : %d", m->max_len); */
/*         int print_until =  m->max_len;//MACRO_MIN(m->min_len, 16); */
/*         for(int i = 0 ; i <  m->seq_comp->L ;i++){ */
/*                 for(int j = 0 ; j < print_until;j++){ */
/*                         fprintf(stdout,"%d ",m->seq_comp->data[i][j]); */
/*                 } */
/*                 fprintf(stdout,"\n"); */
/*         } */

/*         return OK; */
/* ERROR: */
/*         return FAIL; */
/* } */


int collect_seq_comp(struct metrics *m, struct tl_seq *s)
{
        struct seq_composition* c = NULL;
        int mapq_idx = 0;
        int read = 1;
        int start = 0;
        int end = 0;
        RUN(set_read_mapqbin(s, m->mapq_map, &read, &mapq_idx));

        /* s->seq */
        if(read == 1){
                c = m->seq_comp_R1[mapq_idx];
        }else if(read ==2){
                c = m->seq_comp_R2[mapq_idx];
        }
        if(s->data){
                struct aln_data* a = s->data;
                start = a->n_clip5;
                end = a->n_clip3;
        }
        int m_len = MACRO_MIN((int)m->report_max_len, s->len-end);
        int idx = 0;
        for(int i = start;i < m_len;i++){

                char let = s->seq->str[i];
                switch (let) {
                case 'A':
                case 'a':
                        c->data[0][idx]++;
                        break;
                case 'C':
                case 'c':
                        c->data[1][idx]++;
                        break;
                case 'G':
                case 'g':
                        c->data[2][idx]++;
                        break;
                case 'T':
                case 't':
                        c->data[3][idx]++;
                        break;
                default:
                        c->data[4][idx]++;

                        break;
                }
                idx++;
        }
        c->n_counts++;
        return OK;
ERROR:
        return FAIL;
}


int collect_qual_comp(struct metrics *m, struct tl_seq *s, const int offset)
{
        struct qual_composition* c = NULL;
        int mapq_idx = 0;
        int read = 1;
        int start = 0;
        int end = 0;
        RUN(set_read_mapqbin(s, m->mapq_map, &read, &mapq_idx));

        /* s->seq */
        if(read == 1){
                c = m->qual_comp_R1[mapq_idx];
        }else if(read ==2){
                c = m->qual_comp_R2[mapq_idx];
        }
        /* LOG_MSG("Qual len: %d", s->qual->len); */
        /* s->seq */
        if(s->qual->str[0] == 0xFF){

                /* LOG_MSG("Counting : %d", c->n_counts); */
                /* LOG_MSG("Bad qual"); */
                return OK;
        }/* else{ */
        /*         LOG_MSG("Detected a non  missing qual %d", s->qual[0]); */
        /* } */
        if(s->data){
                struct aln_data* a = s->data;
                start = a->n_clip5;
                end = a->n_clip3;
        }

        int m_len = MACRO_MIN((int)m->report_max_len, s->len-end);
        int idx = 0;
        for(int i = start;i < m_len;i++){
                int q = (int)s->qual->str[i] - offset;
                if(q > 41){
                        q = 41;
                }
                c->data[idx][q]++;
                idx++;
        }
        c->n_counts++;

        return OK;
ERROR:
        return FAIL;
}

int collect_error_comp(struct metrics *m, struct tl_seq *s)
{
        if(m->is_aligned == 0){
                return OK;
        }
        struct error_composition* c = NULL;
        struct aln_data* a = NULL;
        uint8_t* g = NULL;
        uint8_t* r = NULL;
        int rp;
        int aln_len;
        int mapq_idx = 0;
        int read = 1;
        RUN(set_read_mapqbin(s, m->mapq_map, &read, &mapq_idx));

        /* s->seq */
        if(read == 1){
                c = m->error_comp_R1[mapq_idx];
        }else if(read ==2){
                c = m->error_comp_R2[mapq_idx];
        }

        a = s->data;
        aln_len = a->aln_len;
        r = a->read;
        g = a->genome;
        /* if(mapq_idx == 3){ */
                /* fprintf(stdout,"%s\n%s\n%s\n",TLD_STR(s->name), r,g); */
        /* } */

        int m_len = MACRO_MIN((int)m->report_max_len, aln_len);
        rp = 0;
        for(int i = 0;i < m_len;i++){
                /* three possibilities */

                if(r[i] == '-' && g[i] == '-'){

                }else if(r[i] != '-' && g[i] == '-'){
                        if(i){
                                if(g[i-1] != '-'){
                                        c->ins[rp]++;
                                        c->n_ins++;
                                }
                        }else{
                                c->ins[rp]++;
                                c->n_ins++;
                        }
                        rp++;
                }else if(r[i] == '-' && g[i] != '-'){
                        if(i){
                                if(r[i-1] != '-'){
                                        c->del[rp]++;
                                        c->n_del++;
                                }
                        }else{
                                c->del[rp]++;
                                c->n_del++;
                        }
                }else if(r[i] != '-' && g[i] != '-'){
                        if(r[i] != g[i]){
                                switch ( r[i]) {
                                case 'A':
                                case 'a':
                                        c->mis[0][rp]++;

                                        break;
                                case 'C':
                                case 'c':
                                        c->mis[1][rp]++;

                                        break;
                                case 'G':
                                case 'g':
                                        c->mis[2][rp]++;
                                        break;
                                case 'T':
                                case 't':
                                        c->mis[3][rp]++;
                                        break;
                                default:
                                        c->mis[4][rp]++;
                                        break;
                                }
                                c->n_mis++;
                        }
                        rp++;

                }
        }
        /* c->n_counts++; */
        return OK;
ERROR:
        return FAIL;
}

int collect_len_comp(struct metrics *m, struct tl_seq *s)
{
        struct len_composition* c = NULL;
        int mapq_idx = 0;
        int read = 1;
        int shorten = 0;
        RUN(set_read_mapqbin(s, m->mapq_map, &read, &mapq_idx));

        if(s->data){
                struct aln_data* a = s->data;
                shorten = a->n_clip5 + a->n_clip3;
        }


        /* s->seq */
        if(read == 1){
                c = m->len_comp_R1[mapq_idx];
        }else if(read ==2){
                c = m->len_comp_R2[mapq_idx];
        }

        c->data[s->len - shorten]++;
        c->n_counts++;
        return OK;
ERROR:
        return FAIL;
}

int set_read_mapqbin( struct tl_seq *s, struct mapqual_bins* map, int *read, int *idx)
{
        *idx = MAPQUALBIN_UNMAP;
        *read = 1;
        if(s->data){
                struct aln_data* a = NULL;
                a = s->data;
                if(a->reverse){
                        LOG_MSG("Oh no I should reverse!");
                        rev_comp_tl_seq(s);
                        if(a->aln_len){
                                reverse_comp(a->genome, a->aln_len);
                                reverse_comp(a->read, a->aln_len);
                        }
                        a->reverse = 0;
                }

                *idx = map->map[ a->map_q];
                if(a->flag  & BAM_FREAD1){
                        *read = 1;
                }else if(a->flag  & BAM_FREAD2){
                        *read = 2;
                }

                if(a->flag & BAM_FUNMAP){
                        *idx = MAPQUALBIN_UNMAP;
                }
        }
        return OK;
}

int metrics_alloc(struct metrics **metrics, int report_max_len)
{

        struct metrics* m = NULL;

        MMALLOC(m, sizeof(struct metrics));
        m->seq_comp_R1 = NULL;
        m->qual_comp_R1 = NULL;
        m->error_comp_R1 = NULL;
        m->len_comp_R1 = NULL;

        m->seq_comp_R2 = NULL;
        m->qual_comp_R2 = NULL;
        m->error_comp_R2 = NULL;
        m->len_comp_R2 = NULL;

        m->mapq_map = NULL;

        m->min_len_R1 = INT_MAX;
        m->max_len_R1 = 0;
        m->min_len_R2 = INT_MAX;
        m->max_len_R2 = 0;
        m->n_R1_reads = 0;
        m->n_R2_reads = 0;

        m->report_max_len = 500;
        m->n_proper_paired = 0;

        if(report_max_len != -1){
                m->report_max_len = report_max_len;
        }

        m->n_paired = 0;

        RUN(get_mapqual_bins(&m->mapq_map));
        m->n_mapq_bins = m->mapq_map->n_bin;
        MMALLOC(m->seq_comp_R1, sizeof(struct seq_composition*) * m->mapq_map->n_bin);
        MMALLOC(m->qual_comp_R1, sizeof(struct qual_composition*) * m->mapq_map->n_bin);
        MMALLOC(m->error_comp_R1, sizeof(struct error_composition*) * m->mapq_map->n_bin);
        MMALLOC(m->len_comp_R1, sizeof(struct len_composition*) * m->mapq_map->n_bin);

        MMALLOC(m->seq_comp_R2, sizeof(struct seq_composition*) * m->mapq_map->n_bin);
        MMALLOC(m->qual_comp_R2, sizeof(struct qual_composition*) * m->mapq_map->n_bin);
        MMALLOC(m->error_comp_R2, sizeof(struct error_composition*) * m->mapq_map->n_bin);
        MMALLOC(m->len_comp_R2, sizeof(struct len_composition*) * m->mapq_map->n_bin);


        for(int i = 0; i < m->mapq_map->n_bin;i++){
                m->seq_comp_R1[i] = NULL;
                m->qual_comp_R1[i] = NULL;
                m->error_comp_R1[i] = NULL;
                m->len_comp_R1[i] = NULL;

                m->seq_comp_R2[i] = NULL;
                m->qual_comp_R2[i] = NULL;
                m->error_comp_R2[i] = NULL;
                m->len_comp_R2[i] = NULL;
        }

        *metrics = m;
        return OK;
ERROR:
        metrics_free(m);
        return FAIL;
}

void metrics_free(struct metrics *m)
{
        if(m){

                for(int i = 0; i < m->n_mapq_bins;i++){
                        if(m->seq_comp_R1[i]){
                                free_seq_comp(m->seq_comp_R1[i]);
                        }
                        if(m->qual_comp_R1[i]){
                                free_qual_comp(m->qual_comp_R1[i]);
                        }
                        if(m->error_comp_R1[i]){
                                free_error_comp(m->error_comp_R1[i]);
                        }
                        if(m->len_comp_R1[i]){
                                free_len_comp(m->len_comp_R1[i]);
                        }

                        if(m->seq_comp_R2[i]){
                                free_seq_comp(m->seq_comp_R2[i]);
                        }
                        if(m->qual_comp_R2[i]){
                                free_qual_comp(m->qual_comp_R2[i]);
                        }
                        if(m->error_comp_R2[i]){
                                free_error_comp(m->error_comp_R2[i]);
                        }
                        if(m->len_comp_R2[i]){
                                free_len_comp(m->len_comp_R2[i]);
                        }

                }
                MFREE(m->seq_comp_R1);
                MFREE(m->qual_comp_R1);
                MFREE(m->error_comp_R1);
                MFREE(m->len_comp_R1);

                MFREE(m->seq_comp_R2);
                MFREE(m->qual_comp_R2);
                MFREE(m->error_comp_R2);
                MFREE(m->len_comp_R2);

                if(m->mapq_map){
                        free_mapqual_bins(m->mapq_map);
                }

                MFREE(m);
        }
}

int alloc_qual_comp(struct qual_composition **qual_comp, int report_max_len)
{
        struct qual_composition* q = NULL;

        MMALLOC(q, sizeof(struct qual_composition));
        q->len = report_max_len;
        q->L = 42;
        q->n_counts = 0;
        q->data = NULL;


        galloc(&q->data, q->len,q->L);
        for(int i = 0;i < q->len; i++){
                for(int j = 0 ; j < q->L;j++){
                        q->data[i][j] = 0;
                }
        }
        *qual_comp = q;
        return OK;
ERROR:
        return FAIL;
}

/* int resize_qual_comp(struct qual_composition* q, int new_max_len) */
/* { */

/*         ASSERT(new_max_len >= q->len,"New len is shorter??"); */

/*         uint32_t** new = NULL; */

/*         galloc(&new, new_max_len,q->L ); */

/*         for(int i = 0; i < new_max_len ;i++){ */
/*                 for(int j = 0; j < q->L;j++){ */
/*                         new[i][j] = 0; */
/*                 } */
/*         } */
/*         for(int i = 0;i < q->len;i++){ */
/*                 for(int j = 0; j < q->L;j++){ */
/*                         new[i][j] = q->data[i][j]; */
/*                 } */
/*         } */
/*         gfree(q->data); */
/*         q->data = new; */
/*         q->len = new_max_len; */
/*         return OK; */
/* ERROR: */
/*         return FAIL; */
/* } */

void free_qual_comp(struct qual_composition *q)
{
        if(q){
                if(q->data){
                        gfree(q->data);
                }
                MFREE(q);
        }
}

int alloc_error_comp(struct error_composition** error_comp, int L,int report_max_len)
{
        struct error_composition* e = NULL;
        MMALLOC(e, sizeof(struct error_composition));
        e->len = report_max_len;
        if(L == TL_SEQ_BUFFER_DNA){
                e->L = 5;
        }else{
                e->L = 21;
        }
        e->n_mis = 0;
        e->n_ins = 0;
        e->n_del = 0;

        e->mis = NULL;
        e->ins = NULL;
        e->del = NULL;

        galloc(&e->mis,e->L, e->len);
        galloc(&e->ins,e->len);
        galloc(&e->del,e->len);

        for(int i = 0; i < e->L;i++){
                for(int j = 0; j < e->len;j++){
                        e->mis[i][j] = 0;
                }
        }

        for(int i = 0; i < e->len;i++){
                e->ins[i] = 0;
                e->del[i] = 0;
        }


        *error_comp = e;
        return OK;
ERROR:
        free_error_comp(e);
        return FAIL;
}

/* int resize_error_comp(struct error_composition *e, int newL, int new_max_len) */
/* { */
/*         /\* ASSERT(newL == e->L,"Alphabet changed???"); *\/ */
/*         ASSERT(new_max_len >= e->len,"New len is shorter??"); */

/*         uint32_t** new = NULL; */
/*         uint32_t* tmp = NULL; */

/*         if(newL == TL_SEQ_BUFFER_DNA){ */
/*                 newL = 5; */
/*         }else{ */
/*                 newL = 21; */
/*         } */
/*         galloc(&new, newL, new_max_len); */

/*         for(int i = 0; i < newL ;i++){ */
/*                 for(int j = 0; j < new_max_len;j++){ */
/*                         new[i][j] = 0; */
/*                 } */
/*         } */
/*         for(int i = 0;i < e->L;i++){ */
/*                 for(int j = 0; j < e->len;j++){ */
/*                         new[i][j] = e->mis[i][j]; */
/*                 } */
/*         } */
/*         gfree(e->mis); */

/*         e->mis = new; */

/*         galloc(&tmp,new_max_len); */
/*         for(int i = 0; i < new_max_len;i++){ */
/*                 tmp[i] = 0; */
/*         } */
/*         for(int i = 0; i < e->len;i++){ */
/*                 tmp[i] = e->ins[i]; */
/*         } */

/*         gfree(e->ins); */
/*         e->ins = tmp; */
/*         tmp = NULL; */

/*         galloc(&tmp,new_max_len); */
/*         for(int i = 0; i < new_max_len;i++){ */
/*                 tmp[i] = 0; */
/*         } */
/*         for(int i = 0; i < e->len;i++){ */
/*                 tmp[i] = e->del[i]; */
/*         } */

/*         gfree(e->del); */
/*         e->del = tmp; */
/*         tmp = NULL; */

/*         e->L  = newL; */
/*         e->len = new_max_len; */
/*         return OK; */
/* ERROR: */
/*         return FAIL; */

/* } */

void free_error_comp(struct error_composition *e)
{
        if(e){
                if(e->mis){
                        gfree(e->mis);
                }
                if(e->ins){
                        gfree(e->ins);
                }
                if(e->del){
                        gfree(e->del);
                }
                MFREE(e);
        }

}


int alloc_seq_comp(struct seq_composition** seq_comp,int L,int report_max_len)
{
        struct seq_composition* s = NULL;
        MMALLOC(s, sizeof(struct seq_composition));

        s->data = NULL;
        s->len = report_max_len;

        if(L == TL_SEQ_BUFFER_DNA){
                s->L = 5;
        }else{
                s->L = 21;
        }

        galloc(&s->data, s->L, s->len);
        for(int i = 0; i <  s->L;i++){
                for(int j = 0; j < s->len;j++){
                        s->data[i][j] = 0;
                }
        }
        s->n_counts = 0;
        *seq_comp = s;

        return OK;
ERROR:
        return FAIL;
}

/* int resize_seq_comp(struct seq_composition* s,int newL, int new_max_len) */
/* { */
/*         /\* ASSERT(newL == s->L,"Alphabet changed???"); *\/ */
/*         /\* ASSERT(new_max_len >= s->len,"New len is shorter??"); *\/ */

/*         /\* uint32_t** new = NULL; *\/ */

/*         /\* if(newL == TL_SEQ_BUFFER_DNA){ *\/ */
/*         /\*         newL = 5; *\/ */
/*         /\* }else{ *\/ */
/*         /\*         newL = 21; *\/ */
/*         /\* } *\/ */


/*         /\* galloc(&new, newL, new_max_len); *\/ */

/*         /\* for(int i = 0; i < newL ;i++){ *\/ */
/*         /\*         for(int j = 0; j < new_max_len;j++){ *\/ */
/*         /\*                 new[i][j] = 0; *\/ */
/*         /\*         } *\/ */
/*         /\* } *\/ */
/*         /\* /\\* LOG_MSG("%d %d -> %d %d", s->L, s->len, newL, new_max_len); *\\/ *\/ */
/*         /\* for(int i = 0;i < s->L;i++){ *\/ */
/*         /\*         for(int j = 0; j < s->len;j++){ *\/ */
/*         /\*                 new[i][j] = s->data[i][j]; *\/ */
/*         /\*         } *\/ */
/*         /\* } *\/ */
/*         /\* gfree(s->data); *\/ */

/*         /\* s->data = new; *\/ */

/*         /\* /\\* s->L  = newL; *\\/ *\/ */
/*         /\* s->len = new_max_len; *\/ */

/*         return OK; */
/* ERROR: */
/*         return FAIL; */
/* } */


void free_seq_comp(struct seq_composition *s)
{
        if(s){
                if(s->data){
                        gfree(s->data);
                }
                MFREE(s);
        }
}

int alloc_len_comp(struct len_composition** len_comp,int cur_max_len)
{
        struct len_composition* s = NULL;
        MMALLOC(s, sizeof(struct len_composition));
        int alloc_len = 128;
        s->data = NULL;

        while(cur_max_len >= alloc_len){
                alloc_len = alloc_len + alloc_len / 2;
        }
        /* LOG_MSG("Allocating :%d", alloc_len); */
        s->len = alloc_len;


        galloc(&s->data, s->len);
        for(int i = 0; i <  s->len;i++){
                s->data[i] = 0;
        }
        s->n_counts = 0;
        *len_comp = s;

        return OK;
ERROR:
        return FAIL;
}

int resize_len_comp(struct len_composition* s,int new_max_len)
{
        /* ASSERT(newL == s->L,"Alphabet changed???"); */
        /* ASSERT(new_max_len >= s->len,"New len is shorter??"); */
        int alloc_len;
        if(new_max_len > s->len){
                uint32_t* new = NULL;
                alloc_len = s->len;
                while(new_max_len >= alloc_len){
                        alloc_len = alloc_len + alloc_len / 2;
                }

                galloc(&new, alloc_len);

                /* for(int i = 0; i < newL ;i++){ */
                for(int i = 0; i < alloc_len;i++){
                        new[i] = 0;
                }
                for(int i = 0; i < s->len;i++){
                        new[i] = s->data[i];
                }

                gfree(s->data);
                s->data = new;
                s->len = alloc_len;
        }
        return OK;
ERROR:
        return FAIL;
}

void free_len_comp(struct len_composition *s)
{
        if(s){
                if(s->data){
                        gfree(s->data);
                }
                MFREE(s);
        }
}
