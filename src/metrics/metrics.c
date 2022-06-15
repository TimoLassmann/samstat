

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


static int collect_seq_comp(struct metrics *m, struct tl_seq *s);
static int collect_qual_comp(struct metrics *m, struct tl_seq *s, const int offset);
static int collect_error_comp(struct metrics *m, struct tl_seq *s);

static int set_read_mapqbin( struct tl_seq *s, struct mapqual_bins* map, int *read, int *idx);

static int alloc_seq_comp(struct seq_composition **seq_comp, int L, int max_len);
static int resize_seq_comp(struct seq_composition* s,int newL, int new_max_len);
static void free_seq_comp(struct seq_composition *s);

static int alloc_qual_comp(struct qual_composition **qual_comp, int max_len);
static int resize_qual_comp(struct qual_composition* q, int new_max_len);
static void free_qual_comp(struct qual_composition *q);

static int alloc_error_comp(struct error_composition** error_comp, int L, int max_len);
static int resize_error_comp(struct error_composition *e, int newL, int new_max_len);
static void free_error_comp(struct error_composition *e);

static int get_mapqual_bins(struct mapqual_bins **map);
static void free_mapqual_bins(struct mapqual_bins *m);


int get_metrics(struct tl_seq_buffer *sb, struct metrics *m)
{
        /*  */
        /* Sanity checks  */

        int len_change = 0;
        for(int i = 0; i < sb->num_seq;i++){

                /* if(m->max_len < (uint32_t)sb->sequences[i]->len){ */
                /*         m->max_len = sb->sequences[i]->len; */
                /*         len_change = 1; */
                /* } */
                /* if(m->min_len > (uint32_t)sb->sequences[i]->len){ */
                /*         m->min_len = sb->sequences[i]->len; */
                /* } */
                if(sb->sequences[i]->data){

                        struct aln_data* a = NULL;
                        a = sb->sequences[i]->data;
                        if(a->flag & BAM_FPAIRED){
                                m->is_paired = 1;
                        }
                        if(a->flag & BAM_FREAD1){
                                if(m->max_len_R1 < (uint32_t)sb->sequences[i]->len){
                                        m->max_len_R1 = sb->sequences[i]->len;
                                        len_change = 1;
                                }
                                if(m->min_len_R1 > (uint32_t)sb->sequences[i]->len){
                                        m->min_len_R1 = sb->sequences[i]->len;
                                }
                        }else if(a->flag & BAM_FREAD2){
                                if(m->max_len_R2 < (uint32_t)sb->sequences[i]->len){
                                        m->max_len_R2 = sb->sequences[i]->len;
                                        len_change = 2;
                                }
                                if(m->min_len_R2 > (uint32_t)sb->sequences[i]->len){
                                        m->min_len_R2 = sb->sequences[i]->len;
                                }
                        }else{
                                if(m->max_len_R1 < (uint32_t)sb->sequences[i]->len){
                                        m->max_len_R1 = sb->sequences[i]->len;
                                        len_change = 1;
                                }
                                if(m->min_len_R1 > (uint32_t)sb->sequences[i]->len){
                                        m->min_len_R1 = sb->sequences[i]->len;
                                }
                        }
                }else{
                        if(m->max_len_R1 < (uint32_t)sb->sequences[i]->len){
                                m->max_len_R1 = sb->sequences[i]->len;
                                len_change = 1;
                        }
                        if(m->min_len_R1 > (uint32_t)sb->sequences[i]->len){
                                m->min_len_R1 = sb->sequences[i]->len;
                        }
                }
        }
        LOG_MSG("Max len: %d %d", m->max_len_R1, m->max_len_R2);
        /* Sequence composition */
        for(int i = 0; i < m->n_mapq_bins;i++){
                if(!m->seq_comp_R1[i]){
                        RUN(alloc_seq_comp(&m->seq_comp_R1[i], sb->L, m->max_len_R1));
                        RUN(alloc_qual_comp(&m->qual_comp_R1[i], m->max_len_R1));
                        RUN(alloc_error_comp(&m->error_comp_R1[i], sb->L,m->max_len_R1));
                }else if(len_change == 1){   /* need to resize as a longer sequence was found  */
                        LOG_MSG("Changing length: %d L %d ", m->max_len_R1, sb->L);
                        RUN(resize_seq_comp(m->seq_comp_R1[i], sb->L, m->max_len_R1));
                        RUN(resize_qual_comp(m->qual_comp_R1[i], m->max_len_R1));
                        RUN(resize_error_comp(m->error_comp_R1[i], sb->L, m->max_len_R1));
                }
        }
        if(m->is_paired){
                for(int i = 0; i < m->n_mapq_bins;i++){
                        if(!m->seq_comp_R2[i]){
                                RUN(alloc_seq_comp(&m->seq_comp_R2[i], sb->L, m->max_len_R2));
                                RUN(alloc_qual_comp(&m->qual_comp_R2[i], m->max_len_R2));
                                RUN(alloc_error_comp(&m->error_comp_R2[i], sb->L,m->max_len_R2));
                        }else if(len_change == 2){   /* need to resize as a longer sequence was found  */
                                RUN(resize_seq_comp(m->seq_comp_R2[i], sb->L, m->max_len_R2));
                                RUN(resize_qual_comp(m->qual_comp_R2[i], m->max_len_R2));
                                RUN(resize_error_comp(m->error_comp_R2[i], sb->L, m->max_len_R2));
                        }
                }
        }

        for(int i = 0; i < sb->num_seq;i++){
                collect_seq_comp(m, sb->sequences[i]);
                collect_qual_comp(m, sb->sequences[i], sb->base_quality_offset);
                collect_error_comp(m, sb->sequences[i]);
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

        snprintf(m->description[MAPQUALBIN_gt30], 256, "MAPQ >= 30");
        snprintf(m->description[MAPQUALBIN_lt30], 256, "MAPQ < 30");
        snprintf(m->description[MAPQUALBIN_lt10], 256, "MAPQ < 10");
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
        RUN(set_read_mapqbin(s, m->mapq_map, &read, &mapq_idx));

        /* s->seq */
        if(read == 1){
                c = m->seq_comp_R1[mapq_idx];
        }else if(read ==2){
                c = m->seq_comp_R2[mapq_idx];
        }

        for(int i = 0;i < s->len;i++){
                switch ( s->seq[i]) {
                case 'A':
                case 'a':
                        c->data[0][i]++;
                        break;
                case 'C':
                case 'c':
                        c->data[1][i]++;

                        break;
                case 'G':
                case 'g':
                        c->data[2][i]++;
                        break;
                case 'T':
                case 't':
                        c->data[3][i]++;
                        break;
                default:
                        c->data[4][i]++;
                        break;
                }
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
        RUN(set_read_mapqbin(s, m->mapq_map, &read, &mapq_idx));

        /* s->seq */
        if(read == 1){
                c = m->qual_comp_R1[mapq_idx];
        }else if(read ==2){
                c = m->qual_comp_R2[mapq_idx];
        }

        /* s->seq */


        for(int i = 0;i < s->len;i++){
                int q = (int)s->qual[i] - offset;
                if(q > 41){
                        q = 41;
                }
                c->data[i][q]++;
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
        rp = 0;
        for(int i = 0;i < aln_len;i++){
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

int set_read_mapqbin( struct tl_seq *s, struct mapqual_bins* map, int *read, int *idx)
{
        *idx = 0;
        *read = 1;
        if(s->data){
                struct aln_data* a = NULL;
                a = s->data;
                if(a->reverse){
                        /* LOG_MSG("Oh no I should reverse!"); */
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

int metrics_alloc(struct metrics **metrics)
{

        struct metrics* m = NULL;

        MMALLOC(m, sizeof(struct metrics));
        m->seq_comp_R1 = NULL;
        m->qual_comp_R1 = NULL;
        m->error_comp_R1 = NULL;

        m->seq_comp_R2 = NULL;
        m->qual_comp_R2 = NULL;
        m->error_comp_R2 = NULL;

        m->mapq_map = NULL;

        m->min_len_R1 = INT_MAX;
        m->max_len_R1 = 0;
        m->min_len_R2 = INT_MAX;
        m->max_len_R2 = 0;

        m->is_paired = 0;
        m->is_aligned = 0;

        RUN(get_mapqual_bins(&m->mapq_map));
        m->n_mapq_bins = m->mapq_map->n_bin;
        MMALLOC(m->seq_comp_R1, sizeof(struct seq_composition*) * m->mapq_map->n_bin);
        MMALLOC(m->qual_comp_R1, sizeof(struct qual_composition*) * m->mapq_map->n_bin);
        MMALLOC(m->error_comp_R1, sizeof(struct error_composition*) * m->mapq_map->n_bin);

        MMALLOC(m->seq_comp_R2, sizeof(struct seq_composition*) * m->mapq_map->n_bin);
        MMALLOC(m->qual_comp_R2, sizeof(struct qual_composition*) * m->mapq_map->n_bin);
        MMALLOC(m->error_comp_R2, sizeof(struct error_composition*) * m->mapq_map->n_bin);


        for(int i = 0; i < m->mapq_map->n_bin;i++){
                m->seq_comp_R1[i] = NULL;
                m->qual_comp_R1[i] = NULL;
                m->error_comp_R1[i] = NULL;

                m->seq_comp_R2[i] = NULL;
                m->qual_comp_R2[i] = NULL;
                m->error_comp_R2[i] = NULL;
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

                        if(m->seq_comp_R2[i]){
                                free_seq_comp(m->seq_comp_R2[i]);
                        }
                        if(m->qual_comp_R2[i]){
                                free_qual_comp(m->qual_comp_R2[i]);
                        }
                        if(m->error_comp_R2[i]){
                                free_error_comp(m->error_comp_R2[i]);
                        }
                }
                MFREE(m->seq_comp_R1);
                MFREE(m->qual_comp_R1);
                MFREE(m->error_comp_R1);

                MFREE(m->seq_comp_R2);
                MFREE(m->qual_comp_R2);
                MFREE(m->error_comp_R2);

                if(m->mapq_map){
                        free_mapqual_bins(m->mapq_map);
                }

                MFREE(m);
        }
}

int alloc_qual_comp(struct qual_composition **qual_comp, int max_len)
{
        struct qual_composition* q = NULL;

        MMALLOC(q, sizeof(struct qual_composition));
        q->len = max_len;
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

int resize_qual_comp(struct qual_composition* q, int new_max_len)
{

        ASSERT(new_max_len >= q->len,"New len is shorter??");

        uint32_t** new = NULL;

        galloc(&new, new_max_len,q->L );

        for(int i = 0; i < new_max_len ;i++){
                for(int j = 0; j < q->L;j++){
                        new[i][j] = 0;
                }
        }
        for(int i = 0;i < q->len;i++){
                for(int j = 0; j < q->L;j++){
                        new[i][j] = q->data[i][j];
                }
        }
        gfree(q->data);
        q->data = new;
        q->len = new_max_len;
        return OK;
ERROR:
        return FAIL;
}

void free_qual_comp(struct qual_composition *q)
{
        if(q){
                if(q->data){
                        gfree(q->data);
                }
                MFREE(q);
        }
}

int alloc_error_comp(struct error_composition** error_comp, int L, int max_len)
{
        struct error_composition* e = NULL;
        MMALLOC(e, sizeof(struct error_composition));
        e->len = max_len;
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

int resize_error_comp(struct error_composition *e, int newL, int new_max_len)
{
        /* ASSERT(newL == e->L,"Alphabet changed???"); */
        ASSERT(new_max_len >= e->len,"New len is shorter??");

        uint32_t** new = NULL;
        uint32_t* tmp = NULL;

        if(newL == TL_SEQ_BUFFER_DNA){
                newL = 5;
        }else{
                newL = 21;
        }
        galloc(&new, newL, new_max_len);

        for(int i = 0; i < newL ;i++){
                for(int j = 0; j < new_max_len;j++){
                        new[i][j] = 0;
                }
        }
        for(int i = 0;i < e->L;i++){
                for(int j = 0; j < e->len;j++){
                        new[i][j] = e->mis[i][j];
                }
        }
        gfree(e->mis);

        e->mis = new;

        galloc(&tmp,new_max_len);
        for(int i = 0; i < new_max_len;i++){
                tmp[i] = 0;
        }
        for(int i = 0; i < e->len;i++){
                tmp[i] = e->ins[i];
        }

        gfree(e->ins);
        e->ins = tmp;
        tmp = NULL;

        galloc(&tmp,new_max_len);
        for(int i = 0; i < new_max_len;i++){
                tmp[i] = 0;
        }
        for(int i = 0; i < e->len;i++){
                tmp[i] = e->del[i];
        }

        gfree(e->del);
        e->del = tmp;
        tmp = NULL;

        e->L  = newL;
        e->len = new_max_len;
        return OK;
ERROR:
        return FAIL;

}

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


int alloc_seq_comp(struct seq_composition** seq_comp,int L, int max_len)
{
        struct seq_composition* s = NULL;
        MMALLOC(s, sizeof(struct seq_composition));

        s->data = NULL;
        s->len = max_len;
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

int resize_seq_comp(struct seq_composition* s,int newL, int new_max_len)
{
        /* ASSERT(newL == s->L,"Alphabet changed???"); */
        ASSERT(new_max_len >= s->len,"New len is shorter??");

        uint32_t** new = NULL;

        if(newL == TL_SEQ_BUFFER_DNA){
                newL = 5;
        }else{
                newL = 21;
        }


        galloc(&new, newL, new_max_len);

        for(int i = 0; i < newL ;i++){
                for(int j = 0; j < new_max_len;j++){
                        new[i][j] = 0;
                }
        }
        LOG_MSG("%d %d -> %d %d", s->L, s->len, newL, new_max_len);
        for(int i = 0;i < s->L;i++){
                for(int j = 0; j < s->len;j++){
                        new[i][j] = s->data[i][j];
                }
        }
        gfree(s->data);

        s->data = new;

        /* s->L  = newL; */
        s->len = new_max_len;

        return OK;
ERROR:
        return FAIL;
}


void free_seq_comp(struct seq_composition *s)
{
        if(s){
                if(s->data){
                        gfree(s->data);
                }
                MFREE(s);
        }

}
