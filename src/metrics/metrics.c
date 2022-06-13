

#include "H5public.h"
#include "alloc/tld-alloc.h"
#include "sam.h"
#include "tld.h"

#include "../htsinterface/htsglue.h"
#include <limits.h>

#define  METRICS_IMPORT
#include "metrics.h"


static int collect_seq_comp(struct metrics *m, struct tl_seq *s);
static int collect_qual_comp(struct metrics *m, struct tl_seq *s, const int offset);

static int set_read_mapqbin( struct tl_seq *s, struct mapqual_bins* map, int *read, int *idx);

static int alloc_seq_comp(struct seq_composition **seq_comp, int L, int max_len);
static int resize_seq_comp(struct seq_composition* s,int newL, int new_max_len);
static void free_seq_comp(struct seq_composition *s);

static int alloc_qual_comp(struct qual_composition **qual_comp, int max_len);
static int resize_qual_comp(struct qual_composition* q, int new_max_len);
static void free_qual_comp(struct qual_composition *q);

static int get_mapqual_bins(struct mapqual_bins **map);
static void free_mapqual_bins(struct mapqual_bins *m);

static int reverse_complement(struct tl_seq *s);

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
        /* Sequence composition */
        for(int i = 0; i < 6;i++){
                if(!m->seq_comp_R1[i]){
                        RUN(alloc_seq_comp(&m->seq_comp_R1[i], sb->L, m->max_len_R1));
                        RUN(alloc_qual_comp(&m->qual_comp_R1[i], m->max_len_R1));
                }else if(len_change == 1){   /* need to resize as a longer sequence was found  */
                        RUN(resize_seq_comp(m->seq_comp_R1[i], sb->L, m->max_len_R1));
                        RUN(resize_qual_comp(m->qual_comp_R2[i], m->max_len_R1));
                }
        }
        if(m->is_paired){
                for(int i = 0; i < 6;i++){
                        if(!m->seq_comp_R2[i]){
                                RUN(alloc_seq_comp(&m->seq_comp_R2[i], sb->L, m->max_len_R2));
                                RUN(alloc_qual_comp(&m->qual_comp_R2[i], m->max_len_R2));
                        }else if(len_change == 2){   /* need to resize as a longer sequence was found  */
                                RUN(resize_seq_comp(m->seq_comp_R2[i], sb->L, m->max_len_R2));
                                RUN(resize_qual_comp(m->qual_comp_R2[i], m->max_len_R2));
                        }
                }
        }

        for(int i = 0; i < sb->num_seq;i++){
                collect_seq_comp(m, sb->sequences[i]);
                collect_qual_comp(m, sb->sequences[i], sb->base_quality_offset);
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
                if(i == 0){
                        m->map[i] = MAPQUALBIN_UNMAP;
                }else if(i < 10){
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

int set_read_mapqbin( struct tl_seq *s, struct mapqual_bins* map, int *read, int *idx)
{
        *idx = 0;
        *read = 1;
        if(s->data){
                struct aln_data* a = NULL;
                a = s->data;
                if(a->reverse){
                        /* LOG_MSG("Oh no I should reverse!"); */
                        reverse_complement(s);
                        a->reverse  = 0;
                }

                *idx = map->map[ a->map_q];
                if(a->flag  & BAM_FREAD1){
                        *read = 1;
                }else if(a->flag  & BAM_FREAD2){
                        *read = 2;
                }

                if(a->flag & BAM_FUNMAP){
                        *idx = 0;
                }
        }
        return OK;
}


int reverse_complement(struct tl_seq *s)
{
        int c;
        c = 0;
        uint8_t* rev = NULL;
        char* revq = NULL;


        MMALLOC(rev, sizeof(uint8_t)* s->malloc_len);
        MMALLOC(revq, sizeof(char)* s->malloc_len);
        for(int i = s->len -1; i >= 0;i--){
                uint8_t l;
                switch ( s->seq[i]) {
                case 'A':
                case 'a':
                        l = 'T';
                        break;
                case 'C':
                case 'c':
                        l = 'G';
                        break;
                case 'G':
                case 'g':
                        l = 'C';
                        break;
                case 'T':
                case 't':
                        l = 'A';
                        break;
                default:
                        l = 'N';
                        break;
                }
                rev[c] = l;
                revq[c] = s->qual[i];
                c++;
        }
        MFREE(s->seq);
        MFREE(s->qual);


        s->seq = rev;
        s->qual = revq;

        return OK;
ERROR:
        return FAIL;
}



int metrics_alloc(struct metrics **metrics)
{

        struct metrics* m = NULL;

        MMALLOC(m, sizeof(struct metrics));
        m->seq_comp_R1 = NULL;
        m->qual_comp_R1 = NULL;

        m->seq_comp_R2 = NULL;
        m->qual_comp_R2 = NULL;

        m->mapq_map = NULL;

        m->min_len_R1 = INT_MAX;
        m->max_len_R1 = 0;
        m->min_len_R2 = INT_MAX;
        m->max_len_R2 = 0;

        m->is_paired = 0;
        m->is_aligned = 0;
        RUN(get_mapqual_bins(&m->mapq_map));

        MMALLOC(m->seq_comp_R1, sizeof(struct seq_composition*) * 6);
        MMALLOC(m->qual_comp_R1, sizeof(struct qual_composition*) * 6);

        MMALLOC(m->seq_comp_R2, sizeof(struct seq_composition*) * 6);
        MMALLOC(m->qual_comp_R2, sizeof(struct qual_composition*) * 6);


        for(int i = 0; i < 6;i++){
                m->seq_comp_R1[i] = NULL;
                m->qual_comp_R1[i] = NULL;

                m->seq_comp_R2[i] = NULL;
                m->qual_comp_R2[i] = NULL;
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
                if(m->mapq_map){
                        free_mapqual_bins(m->mapq_map);
                }
                for(int i = 0; i < 6;i++){
                        if(m->seq_comp_R1[i]){
                                free_seq_comp(m->seq_comp_R1[i]);
                        }
                        if(m->qual_comp_R1[i]){
                                free_qual_comp(m->qual_comp_R1[i]);
                        }
                        if(m->seq_comp_R2[i]){
                                free_seq_comp(m->seq_comp_R2[i]);
                        }
                        if(m->qual_comp_R2[i]){
                                free_qual_comp(m->qual_comp_R2[i]);
                        }
                }
                MFREE(m->seq_comp_R1);
                MFREE(m->qual_comp_R1);

                MFREE(m->seq_comp_R2);
                MFREE(m->qual_comp_R2);

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
        ASSERT(newL == s->L,"Alphabet changed???");
        ASSERT(new_max_len >= s->len,"New len is shorter??");

        uint32_t** new = NULL;

        galloc(&new, newL, new_max_len);

        for(int i = 0; i < newL ;i++){
                for(int j = 0; j < new_max_len;j++){
                        new[i][j] = 0;
                }
        }
        for(int i = 0;i < s->L;i++){
                for(int j = 0; j < s->len;j++){
                        new[i][j] = s->data[i][j];
                }
        }
        gfree(s->data);

        s->data = new;

        s->L  = newL;
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
