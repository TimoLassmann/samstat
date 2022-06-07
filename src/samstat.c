#include "sambamparse/sam_bam_parse.h"
#include "string/str.h"
#include "tld.h"

#include "htslib/sam.h"
#include "htsinterface/htsglue.h"
#include <stdint.h>

/* #define TEST_BAM "/Users/Timo/tmp/77.bam" */
/* #define TEST_BAM "/home/timo/tmp/SYD-40350604.dedup.realigned.recalibrated_UPFB3.bam" */
/* #define TEST_BAM "/home/timo/tmp/test.sam" */
int main(int argc, char *argv[])
{
        fprintf(stdout,"Hello world\n");
        fprintf(stdout,"HTS version:%d\n", HTS_VERSION);
        int* g = NULL;
        galloc(&g,100);

        for(int i = 0; i < 100;i++){
                g[i] = i;
        }

        gfree(g);
        struct sam_bam_file* f_handle = NULL;
        struct tl_seq_buffer* sb = NULL;

        RUN(alloc_tl_seq_buffer(&sb, 20));
        add_aln_data(sb);

        RUN(open_bam(&f_handle, argv[1]));
        while(1){
                RUN(read_bam_chunk(f_handle, sb));
                if(sb->num_seq == 0){
                        break;
                }
                for(int i = 0 ; i < sb->num_seq;i++){
                        struct aln_data* a = (struct aln_data*) sb->sequences[i]->data;
                        fprintf(stdout,"%s\n%s\n", TLD_STR( sb->sequences[i]->name),
                                sb->sequences[i]->seq
                                /* TLD_STR(a->md), */
                                /* a->error, */
                                /* a->reverse */
                                );
                        parse_alignment(sb->sequences[i]);
                        /* get_readable_cigar(a, NULL); */
                }
                break;
/* sb->num_seq = 0; */
                reset_tl_seq_buffer(sb);
        }
        RUN(close_bam(f_handle));


        clear_aln_data(sb);
        free_tl_seq_buffer(sb);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

