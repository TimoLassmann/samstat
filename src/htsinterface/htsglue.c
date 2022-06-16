#include "core/tld-core.h"
#include "sam.h"
#include "tld.h"
#include <htslib/sam.h>

#define HTSGLUE_IMPORT
#include "htsglue.h"

static void free_aln_data_struct(struct aln_data *a);
static int alloc_aln_data(struct aln_data **aln_d);


int add_aln_data(struct tl_seq_buffer *sb)
{
        ASSERT(sb != NULL,"No sb");
        for(int i = 0; i < sb->malloc_num;i++){
                if(!sb->sequences[i]->data){
                        struct aln_data* f = NULL;
                        RUN(alloc_aln_data(&f));
                        sb->sequences[i]->data = f;
                        f = NULL;
                }
        }
        return OK;
ERROR:
        clear_aln_data(sb);
        return FAIL;
}

int clear_aln_data(struct tl_seq_buffer *sb)
{
        ASSERT(sb != NULL,"No sb");
        for(int i = 0; i < sb->malloc_num;i++){
                if(sb->sequences[i]->data){
                        struct aln_data* a = NULL;
                        a = sb->sequences[i]->data;
                        a->md->len = 0;
                }
        }
        return OK;
ERROR:
        return FAIL;
}


int free_aln_data(struct tl_seq_buffer *sb)
{
        ASSERT(sb != NULL,"No sb");
        for(int i = 0; i < sb->malloc_num;i++){
                if(sb->sequences[i]->data){
                        free_aln_data_struct(sb->sequences[i]->data);
                        sb->sequences[i]->data = NULL;
                }
        }
        return OK;
ERROR:
        return FAIL;
}


int read_bam_chunk(struct sam_bam_file *f_handle, struct tl_seq_buffer *sb)
{
        bam1_t *b = f_handle->b;
        bam_hdr_t *h = f_handle->header;
        int r = 0;

        while ((r = sam_read1(f_handle->in, h, b)) >= 0){
                if(!(b->core.flag & BAM_FSECONDARY)){
                        /* if(b->core.qual >= 0){ */
                        struct tl_seq* s = NULL;
                        struct aln_data* a = NULL;
                        uint8_t * seq =    bam_get_seq(b);
                        uint8_t* qual_ptr = bam_get_qual(b);

                        s = sb->sequences[sb->num_seq];
                        a = (struct aln_data* )s->data;

                        a->flag =  b->core.flag;
                        tld_append(s->name, bam_get_qname(b));
                        /* snprintf(s->name,TL_SEQ_MAX_NAME_LEN,"%s",bam_get_qname(b)); */

                        /* snprintf(sb_ptr->name, MAX_SEQ_NAME_LEN,"%s",bam_get_qname(b)); */
                        a->num_hits = 0;
                        a->map_q = b->core.qual;

                        s->len = b->core.l_qseq;

                        /* read in the sequence... */
                        while(s->len+1 >= s->malloc_len){
                                resize_tl_seq(s);
                        }

                        for (int i = 0; i < s->len; ++i){
                                s->seq[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)];
                        }
                        s->seq[s->len ] = 0;

                        /* LOG_MSG("%d", qual_ptr[0]); */
                        if(qual_ptr[0] == 0xFF){
                                s->qual[0] = 0xFF;
                                s->qual[s->len-1] = 0xFF; /* HACK! - this ensures that we have 255
                                                             at the beginning whether or not the sequences
                                                             are on the forward or reverse strand  */
                                /* s->qual[0] = '*'; */
                                /* s->qual[1] = 0; */
                                /* for (int i = 0; i < s->len; ++i){ */
                                        /* s->qual[i] = 'J'; */
                                /* } */
                                /* s->qual[s->len] = 0 ; */
                        }else{
                                for (int i = 0; i < s->len; ++i){
                                        s->qual[i] = qual_ptr[i] + 33;
                                        /* sb_ptr->base_qual[i] = qual_ptr[i] + 33; */
                                }
                                s->qual[s->len] = 0 ;
                        }
                        /* extra stuff  */
                        if(! (BAM_FUNMAP & b->core.flag)){
                                uint8_t * aux =  bam_aux_get(b,"NM");
                                if(aux){
                                        /* LOG_MSG("NM: %d", bam_aux2i(aux)); */
                                        a->error = bam_aux2i(aux);
                                }
                                aux =  bam_aux_get(b,"MD");
                                if(aux){
                                        /* LOG_MSG("NM: %d", bam_aux2i(aux)); */
                                        /* LOG_MSG("%s",bam_aux2Z(aux)); */
                                        tld_append(a->md,bam_aux2Z(aux));
                                }
                        }

                        a->n_cigar = b->core.n_cigar;
                        if(a->n_cigar > a->n_alloc_cigar){
                                a->n_alloc_cigar = MACRO_MAX(a->n_alloc_cigar + a->n_alloc_cigar / 2, a->n_cigar);
                                gfree(a->cigar);
                                a->cigar = NULL;
                                galloc(&a->cigar, a->n_alloc_cigar);
                        }
                        uint32_t* tmp = bam_get_cigar(b);// ((uint32_t*)((b)->data + (b)->core.l_qname))
                        for(int i = 0; i < a->n_cigar;i++){
                                a->cigar[i] = tmp[i];
                        }

                        a->reverse = bam_is_rev(b);
                        if(s->len > sb->max_len){
                                sb->max_len = s->len;
                        }

                        sb->num_seq++;
                        if(sb->num_seq ==  sb->malloc_num){
                                /* exit(0); */
                                return OK;
                        }
                }
        /* } */
        }
        return OK;
ERROR:
        return FAIL;
}

int open_bam(struct sam_bam_file **file, char *filename)
{
        struct sam_bam_file* f_handle = NULL;

        MMALLOC(f_handle, sizeof(struct sam_bam_file));

        f_handle->filename = filename;
        f_handle->b = NULL;
        f_handle->in = NULL;
        f_handle->header = NULL;

        f_handle->in = sam_open(f_handle->filename, "r");
        f_handle->header = sam_hdr_read(f_handle->in);

        f_handle->b = bam_init1();

        *file = f_handle;
        return OK;
ERROR:
        return FAIL;
}

int close_bam(struct sam_bam_file* f_handle)
{
        int status;

        bam_destroy1(f_handle->b);
        bam_hdr_destroy(f_handle->header);


        status = sam_close(f_handle->in);
        if (status < 0) {
                ERROR_MSG("Error closing input: %s.\n",f_handle->filename);
        }
        MFREE(f_handle);
        return OK;
ERROR:
        return FAIL;
}


int alloc_aln_data(struct aln_data **aln_d)
{
        struct aln_data* a = NULL;
        MMALLOC(a, sizeof(struct aln_data));
        a->md = NULL;
        a->cigar = NULL;
        a->genome = NULL;
        a->read = NULL;
        a->aln_len_alloc = 256;
        a->aln_len = 0;
        a->map_q = 0;
        a->num_hits = 0;
        a->n_cigar = 0;
        a->n_alloc_cigar = 16;
        a->reverse = 0;
        RUN(tld_strbuf_alloc(&a->md, 16));
        galloc(&a->read, a->aln_len_alloc);
        galloc(&a->genome, a->aln_len_alloc);
        galloc(&a->cigar, a->n_alloc_cigar);
        *aln_d = a;
        return OK;
ERROR:
        free_aln_data_struct(a);
        return FAIL;
}


void free_aln_data_struct(struct aln_data *a)
{
        if(a){
                if(a->md){
                        tld_strbuf_free(a->md);
                }

                if(a->genome){
                        gfree(a->genome);
                }
                if(a->read){
                        gfree(a->read);
                }

                if(a->cigar){
                        gfree(a->cigar);
                }
                MFREE(a);
        }
}
