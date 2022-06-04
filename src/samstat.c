#include "htsinterface/htsglue.h"
#include "seq/tld-seq.h"
#include "tld.h"

#include <stdint.h>

#define TEST_BAM "/Users/Timo/tmp/77.bam"


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

        RUN(alloc_tl_seq_buffer(&sb, 10));
        add_aln_data(sb);

        RUN(open_bam(&f_handle, TEST_BAM));
        while(1){
                RUN(read_bam_chunk(f_handle, sb));
                if(sb->num_seq == 0){
                        break;
                }
                for(int i = 0 ; i < sb->num_seq;i++){
                        struct aln_data* a = (struct aln_data*) sb->sequences[i]->data;
                        fprintf(stdout,"%s\n%s\n%d\n", sb->sequences[i]->name,
                                sb->sequences[i]->seq,
                                a->error
                                );
                }
                break;
                sb->num_seq = 0;
        }
        RUN(close_bam(f_handle));


        clear_aln_data(sb);
        free_tl_seq_buffer(sb);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

