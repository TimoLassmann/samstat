#include "tld.h"

#include "sam.h"
#include "hts.h"


#define TEST_BAM "~/tmp/SRR528324_strim.fq.all.bam"
struct sam_bam_entry{
	int64_t* start;
	int64_t* stop;
	uint8_t* base_qual;
	char* sequence;
	char* name;
	int qual;
	int len;
	int max_len;
	int num_hits;
	int max_num_hits;
};

struct sam_bam_file{
	samFile *in;
    bam_hdr_t *header;
	bam1_t *b;
	struct seq_info* si;
	char* file_name;
	struct sam_bam_entry** buffer;

	int read_flag;

	int buffer_size;
	int num_read;
	int total_entries_in_file;

	int read_Q_threshold;
	int multimap_Q_threshold;
	int max_num_hits;

	int64_t* cum_chr_len;
	int64_t total_length;
	int read_mode;
};



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
        sam_hdr_t* h = NULL;

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}



struct sam_bam_file* open_SAMBAMfile(char* name,int buffer_size,int max_num_hits, int read_Q_threshold, int multimap_Q_threshold)
{
        struct sam_bam_file* sb_file = NULL;

        int max_len = 120;
        int i;

        MMALLOC(sb_file,sizeof(struct sam_bam_file));
        sb_file->buffer = NULL;
        sb_file->cum_chr_len = NULL;
        sb_file->header = NULL;
        sb_file->in = NULL;
        sb_file->b = NULL;
        sb_file->multimap_Q_threshold = multimap_Q_threshold;
        sb_file->read_Q_threshold = read_Q_threshold;

        sb_file->max_num_hits = max_num_hits;
        sb_file->si = NULL;
        sb_file->total_entries_in_file = 0;
        sb_file->read_flag = BAM_FQCFAIL | BAM_FSECONDARY;
        sb_file->buffer_size = buffer_size;
        sb_file->file_name = strdup(name);
        sb_file->in = sam_open(name, "r");



        sb_file->header = sam_hdr_read(sb_file->in);
        /* bam_cigar_table */
        /* if (sb_file->header->cigar_tab == NULL) { */
        /*         DPRINTF3("need to make header cigartab\n"); */
        /*         MMALLOC(sb_file->header->cigar_tab , sizeof(int8_t) * 128); */

        /*         for (i = 0; i < 128; ++i){ */
        /*                 sb_file->header->cigar_tab[i] = -1; */
        /*         } */
        /*         for (i = 0; BAM_CIGAR_STR[i]; ++i){ */
        /*                 sb_file->header->cigar_tab[(int)BAM_CIGAR_STR[i]] = i; */
        /*         } */
        /* } */

        sb_file->header->ignore_sam_err = 0;

        MMALLOC(sb_file->cum_chr_len, sizeof(int64_t) *sb_file->header->n_targets );
        sb_file->total_length = (int64_t) sb_file->header->target_len[0];
        sb_file->cum_chr_len[0] = 0;
        for(i = 1; i <sb_file->header->n_targets;i++){
                sb_file->cum_chr_len[i]  = sb_file->cum_chr_len[i-1]  + (int64_t) sb_file->header->target_len[i-1];
                sb_file->total_length +=(int64_t) sb_file->header->target_len[i];
        }

        if(buffer_size){
                MMALLOC(sb_file->buffer,sizeof(struct sam_bam_entry*)* buffer_size);
                for(i = 0; i < buffer_size;i++ ){

                        MMALLOC(sb_file->buffer[i],sizeof(struct sam_bam_entry));


                        MMALLOC(sb_file->buffer[i]->sequence, sizeof(char)* max_len);
                        MMALLOC(sb_file->buffer[i]->base_qual,sizeof(uint8_t)* max_len);
                        MMALLOC(sb_file->buffer[i]->start, sizeof(int64_t) * max_num_hits);
                        MMALLOC(sb_file->buffer[i]->stop, sizeof(int64_t) * max_num_hits);
                        MMALLOC(sb_file->buffer[i]->name, sizeof(char) * 256);
                        sb_file->buffer[i]->max_len = max_len;
                        sb_file->buffer[i]->max_num_hits = max_num_hits;
                        sb_file->buffer[i]->num_hits = 0;
                        sb_file->buffer[i]->len = 0;
                        sb_file->buffer[i]->qual = 0;
                }

        }

        /* RUN(make_si_info_in_sam_bam_file(sb_file)); */

        sb_file->b = bam_init1();


        return sb_file;
ERROR:

        if(sb_file){
                if(sb_file->si){
                        /* free_sequence_info(sb_file->si); */
                }
                if(sb_file->b){
                        bam_destroy1(sb_file->b);
                }
                if(sb_file->header){
                        bam_hdr_destroy(sb_file->header);
                }
                if(sb_file->in){
                        i = sam_close(sb_file->in);
                        if(i < 0) WARNING_MSG("Error closing input: %s.\n",sb_file->file_name);
                }

                if(sb_file->buffer){
                        for(i = 0; i < buffer_size;i++ ){
                                if(sb_file->buffer[i]){
                                        MFREE(sb_file->buffer[i]->sequence);
                                        MFREE(sb_file->buffer[i]->start);
                                        MFREE(sb_file->buffer[i]->stop);
                                }
                                MFREE(sb_file->buffer[i]);
                        }
                        MFREE(sb_file->buffer);
                }

                MFREE(sb_file->cum_chr_len);
                MFREE(sb_file->file_name);
                MFREE(sb_file);
        }

        return NULL;
}
