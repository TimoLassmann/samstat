
#include "core/tld-core.h"
#include "seq/tld-seq.h"
#include "tld.h"
#include "sambamparse/sam_bam_parse.h"

#include "htslib/sam.h"
#include "htsinterface/htsglue.h"
#include "param/param.h"
#include <stdint.h>

#define FILE_TYPE_SAMBAM 0
#define FILE_TYPE_FASTAQ 1
#define FILE_TYPE_UNKNOWN 2

int detect_file_type(char *filename, int* type);
int process_sam_bam_file(struct samstat_param* p, int id);
int process_fasta_fastq_file(struct samstat_param* p, int id);
int debug_seq_buffer_print(struct tl_seq_buffer *sb);
/* #define TEST_BAM "/Users/Timo/tmp/77.bam" */
/* #define TEST_BAM "/home/timo/tmp/SYD-40350604.dedup.realigned.recalibrated_UPFB3.bam" */
/* #define TEST_BAM "/home/timo/tmp/test.sam" */
int main(int argc, char *argv[])
{
        struct samstat_param* param = NULL;
        parse_param(argc, argv, &param);

        fprintf(stdout,"Hello world\n");
        fprintf(stdout,"HTS version:%d\n", HTS_VERSION);
        int* g = NULL;
        galloc(&g,100);

        for(int i = 0; i < 100;i++){
                g[i] = i;
        }

        gfree(g);

        for(int i = 0 ; i < param->n_infile;i++){
                int t = -1;
                RUN(detect_file_type(param->infile[i], &t));
                fprintf(stdout,"%s ->type %d\n", param->infile[i],t);
                if(t == FILE_TYPE_FASTAQ){
                        process_fasta_fastq_file(param,i);
                }else if(t == FILE_TYPE_SAMBAM){
                        process_sam_bam_file(param,i);
                }

        }

        param_free(param);
        exit(0);
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
                        char* cig = NULL;
                        get_readable_cigar(a, &cig);
                        fprintf(stdout,"Cigar: %s\n",cig);
                        fprintf(stdout,"   MD: %s\n",TLD_STR(a->md));
                        gfree(cig);

                        fprintf(stdout, "%s (genome)\n",a->genome);
                        for(int i = 0; i < a->aln_len;i++){
                                if(a->genome[i] == a->read[i]){
                                        fprintf(stdout,"|");
                                }else{
                                        fprintf(stdout," ");
                                }
                        }
                        fprintf(stdout,"\n");

                        fprintf(stdout, "%s (read)\n",a->read);

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

int detect_file_type(char *filename, int* type)
{
        if(tld_file_exists(filename)== OK){
                tld_str f = tld_char_to_str(filename);
                tld_str sam_suffix = tld_char_to_str(".sam");
                tld_str bam_suffix = tld_char_to_str(".bam");
                tld_str fasta_suffix = tld_char_to_str(".fasta");
                tld_str fasta_suffix_short = tld_char_to_str(".fa");
                tld_str fastq_suffix = tld_char_to_str(".fastq");
                tld_str fasta_suffix_gz = tld_char_to_str(".fasta.gz");
                tld_str fastq_suffix_gz = tld_char_to_str(".fastq.gz");

                if(tld_suffix_match(f, sam_suffix)){
                        *type = FILE_TYPE_SAMBAM;
                }else if(tld_suffix_match(f, bam_suffix)){
                        *type = FILE_TYPE_SAMBAM;
                }else if(tld_suffix_match(f, fasta_suffix)){
                        *type = FILE_TYPE_FASTAQ;
                }else if(tld_suffix_match(f, fasta_suffix_short)){
                        *type = FILE_TYPE_FASTAQ;
                }else if(tld_suffix_match(f, fastq_suffix)){
                        *type = FILE_TYPE_FASTAQ;
                }else if(tld_suffix_match(f, fasta_suffix_gz)){
                        *type = FILE_TYPE_FASTAQ;
                }else if(tld_suffix_match(f, fastq_suffix_gz)){
                        *type = FILE_TYPE_FASTAQ;
                }else{
                        *type = FILE_TYPE_UNKNOWN;
                }

        }else{
                *type = FILE_TYPE_UNKNOWN;
                ERROR_MSG("File %s not found!!!!", filename);
        }

        return OK;
ERROR:
        return FAIL;
}
int process_sam_bam_file(struct samstat_param* p, int id)
{
        struct sam_bam_file* f_handle = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct alphabet* a = NULL;

        ASSERT(tld_file_exists(p->infile[id]) == OK,"File: %s does not exists",p->infile[id]);

        RUN(alloc_tl_seq_buffer(&sb, p->buffer_size));
        add_aln_data(sb);
        RUN(create_alphabet(&a, NULL,TLALPHABET_DEFAULT_DNA));
        sb->data = a;

        RUN(open_bam(&f_handle, p->infile[id]));
        while(1){

                RUN(read_bam_chunk(f_handle, sb));
                sb->L = TL_SEQ_BUFFER_DNA;
                sb->base_quality_offset = 33; /* Can I get this from sam - don't need to part of the standard */

                if(sb->num_seq == 0){
                        break;
                }
                /* LOG_MSG("L: %d",sb->L); */
                //debug_seq_buffer_print(sb);
                /* break; */
                clear_aln_data(sb);
                reset_tl_seq_buffer(sb);
        }
        RUN(close_bam(f_handle));

        if(sb->data){
                a = sb->data;
                free_alphabet(a);
        }
        free_aln_data(sb);
        free_tl_seq_buffer(sb);
        return OK;
ERROR:
        return FAIL;
}

int process_fasta_fastq_file(struct samstat_param* p, int id)
{
        struct file_handler *f_handle = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct alphabet* a = NULL;
        ASSERT(tld_file_exists(p->infile[id]) == OK,"File: %s does not exists",p->infile[id]);

        RUN(open_fasta_fastq_file(&f_handle, p->infile[id], TLSEQIO_READ));
        RUN(alloc_tl_seq_buffer(&sb, p->buffer_size));

        while(1){

                RUN(read_fasta_fastq_file(f_handle, &sb, p->buffer_size));

                LOG_MSG("%d", sb->base_quality_offset);
                detect_format(sb);
                if(!sb->data){
                        LOG_MSG("Setting alphabet");
                        if(sb->L == TL_SEQ_BUFFER_DNA){
                                RUN(create_alphabet(&a, NULL,TLALPHABET_DEFAULT_DNA));
                        }else if(sb->L == TL_SEQ_BUFFER_PROTEIN){
                                RUN(create_alphabet(&a, NULL,TLALPHABET_DEFAULT_PROTEIN));
                        }else{
                                ERROR_MSG("Bio polymer size %d not recognized", sb->L);
                        }
                        sb->data = a;
                }
                LOG_MSG("%d", sb->base_quality_offset);
                //total_r+= sb->num_seq;
                LOG_MSG("Finished reading chunk: found %d ",sb->num_seq);
                if(sb->num_seq == 0){
                        break;
                }
                /* LOG_MSG("L: %d",sb->L); */
                //debug_seq_buffer_print(sb);
                /* for(int i = 0; i < sb->num_seq;i++){ */
                /*         fprintf(stdout, "%s\n", sb->sequences[i]->seq); */
                /* } */
                reset_tl_seq_buffer(sb);
        }
        if(sb->data){
                a = sb->data;
                free_alphabet(a);
        }
        free_tl_seq_buffer(sb);
        RUN(close_seq_file(&f_handle));
        return OK;
ERROR:
        return  FAIL;
}

int debug_seq_buffer_print(struct tl_seq_buffer *sb)
{
        for(int i = 0 ; i < sb->num_seq;i++){
                struct aln_data* a = (struct aln_data*) sb->sequences[i]->data;
                fprintf(stdout,"%s\n%s\n", TLD_STR( sb->sequences[i]->name),
                        sb->sequences[i]->seq
                        /* TLD_STR(a->md), */
                        /* a->error, */
                        /* a->reverse */
                        );
                if(a){
                        parse_alignment(sb->sequences[i]);
                        /* get_readable_cigar(a, NULL); */
                        char* cig = NULL;
                        get_readable_cigar(a, &cig);
                        fprintf(stdout,"Cigar: %s\n",cig);
                        fprintf(stdout,"   MD: %s\n",TLD_STR(a->md));
                        gfree(cig);

                        fprintf(stdout, "%s (genome)\n",a->genome);
                        for(int i = 0; i < a->aln_len;i++){
                                if(a->genome[i] == a->read[i]){
                                        fprintf(stdout,"|");
                                }else{
                                        fprintf(stdout," ");
                                }
                        }
                        fprintf(stdout,"\n");

                        fprintf(stdout, "%s (read)\n",a->read);
                }
        }
        return OK;
}
