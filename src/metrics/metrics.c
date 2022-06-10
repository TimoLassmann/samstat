

#include "tld.h"

#include "../htsinterface/htsglue.h"
#include <limits.h>

#define  METRICS_IMPORT
#include "metrics.h"


static int collect_seq_comp(struct metrics *m, struct tl_seq *s);
static int alloc_seq_comp(struct seq_composition **seq_comp, int L, int max_len);
static int resize_seq_comp(struct seq_composition* s,int newL, int new_max_len);
static void free_seq_comp(struct seq_composition *s);

static int get_mapqual_bins(struct mapqual_bins **map);
static void free_mapqual_bins(struct mapqual_bins *m);

static int reverse_complement(struct tl_seq *s);

int get_mapqual_bins(struct mapqual_bins **map)
{
        struct mapqual_bins* m = NULL;

        MMALLOC(m, sizeof(struct mapqual_bins));

        m->map = NULL;
        m->description = NULL;
        m->len = 256;

        galloc(&m->map,m->len);
        galloc(&m->description, 6, 256);

        for(int i = 0; i < m->len;i++){
                if(i == 0){
                        m->map[i] = 5;
                }else if(i < 3){
                        m->map[i] = 4;
                }else if(i < 10){
                        m->map[i] = 3;
                }else if(i < 20){
                        m->map[i] = 2;
                }else if(i < 30){
                        m->map[i] = 1;
                }else{
                        m->map[i] = 0;
                }
        }

        snprintf(m->description[0], 256, "MAPQ >= 30");
        snprintf(m->description[1], 256, "MAPQ < 30");
        snprintf(m->description[2], 256, "MAPQ < 20");
        snprintf(m->description[3], 256, "MAPQ < 10");
        snprintf(m->description[4], 256, "MAPQ < 3");
        snprintf(m->description[5], 256, "Unmapped");

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

int get_metrics(struct tl_seq_buffer *sb, struct metrics *m)
{
        /*  */
        /* Sanity checks  */

        int len_change = 0;
        for(int i = 0; i < sb->num_seq;i++){
                if(m->max_len < (uint32_t)sb->sequences[i]->len){
                        m->max_len = sb->sequences[i]->len;
                        len_change = 1;
                }
                if(m->min_len > (uint32_t)sb->sequences[i]->len){
                        m->min_len = sb->sequences[i]->len;
                }
        }
        /* Sequence composition */
        for(int i = 0; i < 6;i++){
                if(!m->seq_comp[i]){
                        RUN(alloc_seq_comp(&m->seq_comp[i], sb->L, m->max_len));
                }else if(len_change){   /* need to resize as a longer sequence was found  */
                        RUN(resize_seq_comp(m->seq_comp[i], sb->L, m->max_len));
                }
        }

        for(int i = 0; i < sb->num_seq;i++){
                collect_seq_comp(m, sb->sequences[i]);
        }


        return OK;
ERROR:
        return FAIL;
}

int collect_seq_comp(struct metrics *m, struct tl_seq *s)
{
        struct seq_composition* c = NULL;
        int mapq_idx = 0;
        if(s->data){
                struct aln_data* a = NULL;
                a = s->data;
                if(a->reverse){
                        /* LOG_MSG("Oh no I should reverse!"); */
                        reverse_complement(s);
                }
                mapq_idx = m->mapq_map->map[ a->map_q];
        }
        /* s->seq */
        c = m->seq_comp[mapq_idx];

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
        m->seq_comp = NULL;
        m->mapq_map = NULL;
        m->min_len = INT_MAX;
        m->max_len = 0;

        RUN(get_mapqual_bins(&m->mapq_map));

        MMALLOC(m->seq_comp, sizeof(struct seq_composition*) * 6);

        for(int i = 0; i < 6;i++){
                m->seq_comp[i] = NULL;
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
                        if(m->seq_comp[i]){
                                free_seq_comp(m->seq_comp[i]);
                        }

                }
                MFREE(m->seq_comp);
                MFREE(m);
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
