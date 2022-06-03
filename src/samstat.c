#include "htsinterface/htsglue.h"
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
        RUN(open_bam(&f_handle, TEST_BAM));
        RUN(read_bam_chunk(f_handle, NULL));
        RUN(close_bam(f_handle));

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

