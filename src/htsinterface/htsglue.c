#include "seq/tld-seq.h"
#include "tld.h"

#define HTSGLUE_IMPORT
#include "htsglue.h"

int read_bam_chunk(struct sam_bam_file *f_handle, struct tl_seq_buffer *sb)
{
        bam1_t *b = f_handle->b;
        bam_hdr_t *h = f_handle->header;
        int r = 0;

        while ((r = sam_read1(f_handle->in, h, b)) >= 0){
                if(b->core.qual >= 1){
                        uint8_t * seq =    bam_get_seq(b);
                        uint8_t* qual_ptr = bam_get_qual(b);
                        int len =b->core.l_qseq;
                        for (int i = 0; i < len; ++i){
                                fprintf(stdout,"%c","=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)]);
                        }
                        fprintf(stdout,"\n");
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
