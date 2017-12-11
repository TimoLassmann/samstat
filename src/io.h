#define SEEK_START 0
#define SEEK_END 2

#include <unistd.h>

struct read_info{
	char* name;
	char* qual;
	char* seq;
	char* labels;
	char* cigar;
	char* md;
	int errors;
	float mapq;
	int len;
	int strand;
};


FILE* io_handler(FILE* file, int file_num,struct parameters* param);
void print_seq(struct read_info* ri,FILE* out);
int read_sam_chunk(struct read_info** ri,struct parameters* param,FILE* file);
int read_fasta_fastq(struct read_info** ri,struct parameters* param,FILE *file);

int qsort_ri_mapq_compare(const void *a, const void *b);


int file_exists (char * name);


struct read_info** malloc_read_info(struct read_info** ri, int numseq);
struct read_info** clear_read_info(struct read_info** ri, int numseq);
void free_read_info(struct read_info** ri, int numseq);



