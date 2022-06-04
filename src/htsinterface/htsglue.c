#include "alloc/tld-alloc.h"
#include "core/tld-core.h"
#include "seq/tld-seq.h"
#include "tld.h"
#include <stdio.h>

#define HTSGLUE_IMPORT
#include "htsglue.h"

static void free_aln_data(struct aln_data *a);
static int alloc_aln_data(struct aln_data **aln_d);


int add_aln_data(struct tl_seq_buffer *sb)
{
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
        for(int i = 0; i < sb->malloc_num;i++){
                if(sb->sequences[i]->data){
                        free_aln_data(sb->sequences[i]->data);
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
                if(b->core.qual >= 1){
                        struct tl_seq* s = NULL;
                        struct aln_data* a = NULL;
                        uint8_t * seq =    bam_get_seq(b);
                        uint8_t* qual_ptr = bam_get_qual(b);
                        int len =b->core.l_qseq;

                        s = sb->sequences[sb->num_seq];
                        a = (struct aln_data* )s->data;


                        for (int i = 0; i < len; ++i){
                                fprintf(stdout,"%c","=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)]);
                        }
                        fprintf(stdout,"\n");
                        snprintf(s->name,TL_SEQ_MAX_NAME_LEN,"%s",bam_get_qname(b));

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

                        if(qual_ptr[0] == 0xFF){
                                s->qual[0] = '*';
                                s->qual[1] = 0;
                        }else{
                                for (int i = 0; i < s->len; ++i){
                                        s->qual[i] = qual_ptr[i] + 33;
                                        /* sb_ptr->base_qual[i] = qual_ptr[i] + 33; */
                                }
                                s->qual[s->len] = 0 ;
                        }
                        /* extra stuff  */
                        if(! (BAM_FUNMAP & b->core.flag)){
                                uint8_t * aux =  bam_aux_get(b,"NH");
                                LOG_MSG("NM: %d", bam_aux2i(aux));
                                a->error = bam_aux2i(aux);

                        }

                        sb->num_seq++;
                        if(sb->num_seq ==  sb->malloc_num){
                                /* exit(0); */
                                return OK;
                        }
                }
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
        int i;
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
        a->stop = 0;
        a->start = 0;
        a->map_q = 0;
        a->num_hits = 0;

        *aln_d = a;
        return OK;
ERROR:
        free_aln_data(a);
        return FAIL;
}


void free_aln_data(struct aln_data *a)
{
        if(a){
                MFREE(a);
        }
}
