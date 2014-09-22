/*
 
 Copyright (C) 2013 Timo Lassmann <timolassmann@gmail.com>
 
 This file is part of TagDust.
 
 TagDust is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 TagDust is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Tagdust.  If not, see <http://www.gnu.org/licenses/>.
 
 */

#ifndef MMALLOC
#include "malloc_macro.h"
#endif

#include "interface.h"
#include "nuc_code.h"
#include "misc.h"
#include "io.h"
#include "hmm.h"
#include "viz.h"
#include "pst.h"

#include <ctype.h>
#include <time.h>
#include <sys/stat.h>


#define MAX_SEQ_LEN 1000
#define MAXERROR 100

struct seq_stats{
	int** seq_len;
	int*** nuc_composition;
	int** seq_quality;
	int** seq_quality_count;
	int* aln_quality;
	int* alignments;
	int* nuc_num;
	float** errors;
	double* percent_identity;
	int*** mismatches;
	int*** insertions;
	int** deletions;
	int* base_qualities;
	int base_quality_offset;
	int sam;
	int md;
	int min_len;
	int max_len;
	int average_len;
	int max_error_per_read;
	int total_reads;
};


struct hmm* init_samstat_hmm(int average_length, int max_sequence_len);
struct seq_stats* init_seq_stats(void);
struct seq_stats* clear_seq_stats(struct seq_stats* seq_stats);
struct seq_stats* reformat_base_qualities(struct seq_stats* seq_stats);

void free_seq_stats(struct seq_stats* seq_stats);
void print_stats(struct seq_stats* seq_stats);
int parse_cigar_md(struct read_info* ri,struct seq_stats* seq_stats,int qual_key);

char* make_file_stats(char* filename,char* buffer);

int main (int argc,char * argv[])
{
	struct parameters* param = NULL;
	struct seq_stats* seq_stats = NULL;
	struct plot_data* pd = NULL;
	struct hmm_data* hmm_data= NULL;
	struct hmm** hmms = NULL;
	
	int (*fp)(struct read_info** ,struct parameters*,FILE* ) = NULL;
	FILE* file = NULL;
	FILE* outfile = NULL;
	int numseq = 0;
	int i,j,c,fileID;
	int qual_key = 0;
	int aln_len = 0;
	int first_lot =1;

	
	init_nuc_code();
	
	param = interface(param,argc,argv);
	
	
#ifdef DEBUG
	param->num_query = 100;
#else
	param->num_query = 1000000;
#endif
	struct read_info** ri = NULL;
	ri = malloc_read_info(ri, param->num_query );
	
	seq_stats = init_seq_stats();
	
	MMALLOC(hmm_data, sizeof(struct hmm_data));
	hmm_data->length = 0;
	hmm_data->score = 0;
	hmm_data->string = 0;
	hmm_data->iterations = 5;
	hmm_data->run_mode = MODE_BAUM_WELCH;
	hmm_data->num_threads = 4;
	hmm_data->weight = 0;
	
	MMALLOC(hmm_data->length,sizeof(int) *param->num_query);
	MMALLOC(hmm_data->weight,sizeof(float) *param->num_query);
	MMALLOC(hmm_data->score,sizeof(float) *param->num_query);
	MMALLOC(hmm_data->string , sizeof(char* ) * param->num_query);
	for(i = 0; i < param->num_query;i++){
		hmm_data->length[i] = 0;
		hmm_data->score[i] = prob2scaledprob(0.0);
		hmm_data->weight[i] = prob2scaledprob(1.0);
		hmm_data->string[i] = 0;
	}
	
	for(fileID = 0; fileID < param->infiles;fileID++){
		sprintf(param->buffer,"%s\n--------------------------------------------------\n", shorten_pathname(param->infile[fileID]));
		param->messages = append_message(param->messages, param->buffer);
		seq_stats =  clear_seq_stats(seq_stats);
		sprintf(param->outfile,"%s.samstat.html", param->infile[fileID]);
		
		if ((outfile = fopen(param->outfile, "w")) == NULL){
			sprintf(param->buffer,"ERROR: Cannot open outpuf file: %s\n",param->outfile);
			param->messages = append_message(param->messages, param->buffer);
			//fprintf(stderr,"can't open output\n");
			free_param(param);
			exit(EXIT_FAILURE);
		}

		//outfile
		
		file =  io_handler(file, fileID,param);
		if(param->sam == 0){
			fp = &read_fasta_fastq;
		}else {
			fp = &read_sam_chunk;
		}
		
		seq_stats->sam = param->sam;
		qual_key = 0;
		aln_len = 0;
		first_lot =1;
		
		sprintf(param->buffer,"Starting to collect data.\n");
		param->messages = append_message(param->messages, param->buffer);
		
		while ((numseq = fp(ri, param,file)) != 0){
			for(i = 0; i < numseq;i++){
				if(ri[i]->len > seq_stats->max_len){
					seq_stats->max_len = ri[i]->len;
				}
				if(ri[i]->len < seq_stats->min_len ){
					seq_stats->min_len = ri[i]->len;
				}
				seq_stats->average_len += ri[i]->len;
				
				if(ri[i]->mapq == 0){
					qual_key = 5;
				}else if(ri[i]->mapq < 3){
					qual_key = 0;
				}else if(ri[i]->mapq < 10){
					qual_key = 1;
				}else if(ri[i]->mapq < 20){
					qual_key = 2;
				}else if(ri[i]->mapq < 30){
					qual_key = 3;
				}else{
					qual_key = 4;
				}
				if(ri[i]->cigar && ri[i]->md){
					if(ri[i]->cigar[0] != '*'){
						aln_len = parse_cigar_md(ri[i],seq_stats, qual_key);
						seq_stats->md = 1;
					}
				}
				if(ri[i]->strand != 0){
					ri[i]->seq = reverse_complement(ri[i]->seq,ri[i]->len);
					if(ri[i]->qual[0] != '*'){
						reverse_sequence(ri[i]->qual, ri[i]->len);
					}
				}
				if(ri[i]->qual){
					if(ri[i]->qual[0] != '*'){
						for(j = 0; j < ri[i]->len;j++){
							seq_stats->seq_quality[qual_key][j] += (int)(ri[i]->qual[j]);
							seq_stats->seq_quality_count[qual_key][j] += 1; //(int)(ri[i]->qual[j]);
							seq_stats->base_qualities[(int)(ri[i]->qual[j])]++;
						}
					}
				}
				seq_stats->alignments[qual_key]++;
				seq_stats->total_reads++;
				// sequence length
				if(ri[i]->len >=  MAX_SEQ_LEN){
					seq_stats->seq_len[qual_key][MAX_SEQ_LEN-1]++;
				}else{
					seq_stats->seq_len[qual_key][ri[i]->len]++;
				}
				// sequence composition
				for(j = 0;j <  ri[i]->len;j++){
					seq_stats->nuc_num[(int)ri[i]->seq[j]]++;
					seq_stats->nuc_composition[qual_key][j][(int)ri[i]->seq[j]]++;
				}
				
				if(ri[i]->errors != -1){
					if(ri[i]->errors > seq_stats->max_error_per_read){
						seq_stats->max_error_per_read = ri[i]->errors;
					}
					seq_stats->percent_identity[qual_key] +=(((double)aln_len - (double)ri[i]->errors) / (double)aln_len * 100.0);
					if(ri[i]->errors >= MAXERROR){
						seq_stats->errors[qual_key][MAXERROR-1]++;
					}else{
						seq_stats->errors[qual_key][ri[i]->errors]++;
					}
				}
			}
			//needs to be run after sequences are reverse complemented.....
			if(first_lot){
				first_lot = 0;
				seq_stats->average_len = (int) floor((double) seq_stats->average_len / (double) numseq   + 0.5);
				seq_stats = reformat_base_qualities(seq_stats);
				
				hmms = NULL;
				MMALLOC(hmms,sizeof(struct hmm*) * 3) ;
				
				for(i =0 ; i < 3;i++){
					hmms[i] = NULL;
					hmms[i] = init_samstat_hmm(seq_stats->min_len, seq_stats->max_len);
				}
				
				// run for Q20-40  and unmapped.
				j = 0;
				for(i = 0; i < numseq;i++){
					if(ri[i]->mapq >= 20){
						hmm_data->length[j] = ri[i]->len;
						hmm_data->string[j] = ri[i]->seq;
						hmm_data->weight[j] = prob2scaledprob(1.0);
						j++;
					}
				}
				hmm_data->num_seq =j;
				if(j > 100){
					sprintf(param->buffer,"Training a HMM on mapq > 20 reads.\n");
					param->messages = append_message(param->messages, param->buffer);
					hmms[0] = run_EM_iterations(hmms[0],hmm_data);
					sprintf(param->buffer,"Done.\n");
					param->messages = append_message(param->messages, param->buffer);
					//print_hmm_parameters(hmms[0]);
					//exit(0);
				}else{
					free_hmm(hmms[0]);
					hmms[0] = 0;
				}
				j = 0;
				for(i = 0; i < numseq;i++){
					if(ri[i]->mapq > 0 && ri[i]->mapq < 20){
						hmm_data->length[j] = ri[i]->len;
						hmm_data->string[j] = ri[i]->seq;
						hmm_data->weight[j] = prob2scaledprob(1.0);
						j++;
					}
				}
				hmm_data->num_seq =j;
				if(j > 100){
					sprintf(param->buffer,"Training a HMM on 0  <= mapq  < 20 reads.\n");
					param->messages = append_message(param->messages, param->buffer);
					hmms[1] = run_EM_iterations(hmms[1],hmm_data);
					sprintf(param->buffer,"Done.\n");
					param->messages = append_message(param->messages, param->buffer);
				}else{
					free_hmm(hmms[1]);
					hmms[1] = 0;
				}
				
				j = 0;
				for(i = 0; i < numseq;i++){
					if(ri[i]->mapq ==  0){
						hmm_data->length[j] = ri[i]->len;
						hmm_data->string[j] = ri[i]->seq;
						hmm_data->weight[j] = prob2scaledprob(1.0);
						j++;
					}
				}
				hmm_data->num_seq =j;
				if(j > 100){
					sprintf(param->buffer,"Training a HMM on unmapped reads.\n");
					param->messages = append_message(param->messages, param->buffer);
					hmms[2] = run_EM_iterations(hmms[2],hmm_data);
					sprintf(param->buffer,"Done.\n");
					param->messages = append_message(param->messages, param->buffer);
				}else{
					free_hmm(hmms[2]);
					hmms[2] = 0;
				}
							
				
				
				/*struct pst* pst = NULL;
				 char alphabet[] = "ACGTN";
				 
				 fprintf(stderr,"PST building... %d\n",numseq);
				 pst = alloc_pst(numseq) ;
				 for(i = 0; i < numseq;i++){
				 for(j = 0; j < ri[i]->len;j++){
				 ri[i]->seq[j] = alphabet[(int)ri[i]->seq[j]];
				 }
				 pst->total_len += ri[i]->len;
				 }
				 
				 cStartClock = clock();
				 
				 if(pst->current_suffix_size < pst->total_len){
				 pst->suffix_array = realloc(pst->suffix_array , sizeof(char*)* (pst->total_len+64));
				 pst->current_suffix_size =  (pst->total_len+64);
				 }
				 
				 c = 0;
				 pst->mean_length = 0.0;
				 for(i = 0; i < numseq;i++){
				 for(j = 0; j < ri[i]->len;j++){//ri[i]->len;j++){
				 pst->suffix_array[c] = ri[i]->seq +j;
				 c++;
				 }
				 pst->mean_length +=  ri[i]->len;
				 }
				 
				 pst->mean_length /= (float)numseq;
				 pst->suffix_len = c;
				 pst->numseq = numseq;
				 pst->rank_array = malloc(sizeof(struct ranks*) * (int)(numseq));
				 
				 for(i = 0; i  < numseq;i++){
				 pst->rank_array[i] = malloc(sizeof(struct ranks));
				 }
				 
				 fprintf(stderr,"%d\t%d\t%d	%f\n",c,pst->total_len, pst->suffix_len*(4 + 4*5),pst->mean_length);
				 qsort(pst->suffix_array, pst->suffix_len, sizeof(char *), qsort_string_cmp);
				 fprintf(stderr,"built SA in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);
				 
				 //init root - removes if statement in recursion...
				 double sum = 0.0;
				 char tmp[MAX_PST_LEN+5];
				 for(i = 0;i < 5;i++){
				 tmp[0] = alphabet[i];
				 tmp[1] = 0;//alphabet[i];
				 c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,1);
				 pst->pst_root->nuc_probability[i] = c;
				 pst->ppt_root->nuc_probability[i] = c;
				 sum+= c;
				 }
				 
				 for(i = 0;i < 5;i++){
				 pst->pst_root->nuc_probability[i] =  pst->pst_root->nuc_probability[i]/ sum;
				 pst->ppt_root->nuc_probability[i] =  pst->ppt_root->nuc_probability[i]/ sum;
				 }
				 
				 pst->pst_root = build_pst(pst,pst->pst_root );
				 pst->ppt_root = build_ppt(pst,pst->ppt_root );
				 
				 fprintf(stderr,"built PST in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);
				 cStartClock = clock();
				 pst->pst_root  = alloc_bit_occ_pst(pst->pst_root , numseq);
				 pst->ppt_root = alloc_bit_occ_pst(pst->ppt_root, numseq);
				 ri =  scan_read_with_pst( ri, pst);
				 fprintf(stderr,"scanned  in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);
				 
				 int num_patterns = 0;
				 num_patterns = count_patterns(pst->pst_root, num_patterns);
				 num_patterns = count_patterns(pst->ppt_root,num_patterns);
				 
				 struct pst_node** all_patterns = 0;
				 
				 all_patterns = malloc(sizeof(struct pst_node*)  *num_patterns );
				 
				 num_patterns = 0;
				 num_patterns = add_patterns(all_patterns,pst->pst_root, num_patterns);
				 num_patterns = add_patterns(all_patterns,pst->ppt_root,num_patterns);
				 
				 qsort((void *)  all_patterns, num_patterns, sizeof(struct pst_node* ),(compfn) sort_pst_nodel_according_to_label);
				 
				 fprintf(stderr,"%d numpatterns\n",num_patterns );
				 c = 0;
				 for(i = 0; i < num_patterns-1;i++){
				 if(strcmp(all_patterns[i]->label, all_patterns[i+1]->label)){
				 all_patterns[c] = all_patterns[i];
				 c++;
				 }
				 }
				 
				 qsort((void *)  all_patterns, num_patterns, sizeof(struct pst_node* ),(compfn) sort_pst_nodel_according_to_occ);
				 
				 
				 for(i = 0; i < (num_patterns < 10 ? num_patterns : 10 );i++){
				 fprintf(stderr,"%s	%d\n",all_patterns[i]->label,all_patterns[i]->occ );
				 }*/
			}
		}
		pclose(file);
		
#ifdef DEBUG
		print_stats(seq_stats);
#endif
		
		
		pd = malloc_plot_data(10,   seq_stats->max_len+2   );
		pd->height = 250;
		sprintf(pd->plot_title, "%s",shorten_pathname(param->infile[fileID]));
		pd->description = make_file_stats(param->infile[fileID],pd->description);
		sprintf(param->buffer,"%d reads, %s", seq_stats->total_reads,pd->description);
		sprintf(pd->description,"%s",param->buffer);
		
		
		print_html5_header(outfile,pd);
		
		sprintf(pd->labels[0], "%s","Number");
		sprintf(pd->labels[1], "%s","Percentage");
		
		pd->data[5][0] =  seq_stats->alignments[5];
		pd->data[5][1] =  (float)seq_stats->alignments[5] / (float)seq_stats->total_reads * 100.0;
		sprintf(pd->series_labels[5], "Unmapped");
		
		pd->data[4][0] =  seq_stats->alignments[0];
		pd->data[4][1] =  (float)seq_stats->alignments[0] / (float)seq_stats->total_reads* 100.0;
		sprintf(pd->series_labels[4], "MAPQ  <  3");
		
		pd->data[3][0] =  seq_stats->alignments[1];
		pd->data[3][1] =  (float)seq_stats->alignments[1] / (float)seq_stats->total_reads* 100.0;
		sprintf(pd->series_labels[3], "MAPQ  < 10");
		
		pd->data[2][0] =  seq_stats->alignments[2];
		pd->data[2][1] =  (float)seq_stats->alignments[2] / (float)seq_stats->total_reads* 100.0;
		sprintf(pd->series_labels[2], "MAPQ  < 20");
		
		pd->data[1][0] =  seq_stats->alignments[3];
		pd->data[1][1] =  (float)seq_stats->alignments[3] / (float)seq_stats->total_reads* 100.0;
		sprintf(pd->series_labels[1], "MAPQ  < 30");
		
		pd->data[0][0] =  seq_stats->alignments[4];
		pd->data[0][1] =  (float)seq_stats->alignments[4] / (float)seq_stats->total_reads* 100.0;
		
		sprintf(pd->series_labels[0], "MAPQ >= 30");
		
		pd->num_points = 1;
		pd->num_series = 6;
		pd->color_scheme = 4;
		pd->width = 400;
		//sprintf(pd->description,"");
		sprintf(pd->description,"Number of alignments in various mapping quality (MAPQ) intervals and number of unmapped sequences.");
		
		//pd->description = 0;
		sprintf(pd->plot_title, "Mapping stats:");
		pd->plot_type = PIE_PLOT;
		print_html5_chart(outfile, pd);
		
		sprintf(pd->series_labels[6], "Total");
		pd->data[6][0] =  seq_stats->total_reads;
		pd->data[6][1] =  100.0;
		pd->num_points = 2;
		pd->num_series = 7;
		sprintf(pd->description,"Number of alignments in various mapping quality (MAPQ) intervals and number of unmapped sequences.");
		
		print_html_table(outfile, pd);
		
		
		
		
		pd->color_scheme = 0;
		pd->width = 900;
		
		for(i = 0; i < 6;i++){
			if(!seq_stats->alignments[i]){
				pd->show_series[i] =0;
			}
		}
		
		
		int plots =0;
		
		if(seq_stats->max_len == seq_stats->min_len){
			fprintf(outfile,"<h2>Read Length: All reads are %dnt long</h2>\n",seq_stats->max_len);
			
		}else{
			plots =0;
			for(i = 0; i < 6;i++){
				if(plots ==0){
					sprintf(pd->labels[0], "%dnt",seq_stats->min_len-1);
				}
				pd->data[i][0] = 0 ;
				for(j = seq_stats->min_len; j <= seq_stats->max_len;j++){
					if(plots ==0){
						sprintf(pd->labels[j-seq_stats->min_len+1], "%dnt",j);
					}
					pd->data[i][j-seq_stats->min_len+1] = seq_stats->seq_len[i][j];
				}
			}
			
			pd->width = 700;
			pd->color_scheme = 4;
			pd->num_points = seq_stats->max_len - seq_stats->min_len+2;
			if(pd->num_points < 20){
				pd->num_points_shown =pd->num_points;
			}else{
				pd->num_points_shown = 20;
			}
			pd->num_series = 6;
			sprintf(pd->description,"Distribution of read lengths separated by mapping quality thresholds.");
			sprintf(pd->plot_title, "Read Length Distributions");
			pd->plot_type = LINE_PLOT;
			print_html5_chart(outfile, pd);
			
			//pd->num_points = 0;
			//pd->num_series = 6;
			//sprintf(pd->description,"Distribution of read lengths separated by mapping quality thresholds.");
			
			//print_html_table(stdout, pd);
		}
		// BASE QUALITIES...
		
		if (seq_stats->base_quality_offset != -1){
			for(i = 0; i < 6;i++){
				for(j = 0; j <= seq_stats->max_len;j++){
					if(plots ==0){
						sprintf(pd->labels[j], "%dnt",j+1);
					}
					if(seq_stats->alignments[i] ){
						pd->data[i][j] =  (float)seq_stats->seq_quality[i][j] /   (float)seq_stats->seq_quality_count[i][j] - (float)seq_stats->base_quality_offset;
					}else{
						pd->data[i][j] = 0;
					}
				}
			}
			
			
			pd->width = 700;
			pd->color_scheme = 4;
			pd->num_points = seq_stats->max_len;
			if(pd->num_points < 20){
				pd->num_points_shown =pd->num_points;
			}else{
				pd->num_points_shown = 20;
			}
			pd->num_series = 6;
			sprintf(pd->description,"Base quality distributions separated by mapping quality thresholds.");
			sprintf(pd->plot_title, "Base Quality Distributions");
			pd->plot_type = LINE_PLOT;
			print_html5_chart(outfile, pd);
			
			//pd->num_points = 0;
			//pd->num_series = 6;
			//sprintf(pd->description,"Base quality distributions separated by mapping quality thresholds.");
			
			//print_html_table(stdout, pd);
		}
		
		
		///HMM plots - need to count  nucleotide frequencies.. .
		
		float sum = 0;
		
		for(i = 0; i < 6;i++){
			sum +=seq_stats->nuc_num[i];
		}
		if(hmms[0]){
			pd->color_scheme = 0;
			
			sprintf(pd->series_labels[0],"A");
			sprintf(pd->series_labels[1],"C");
			sprintf(pd->series_labels[2],"G");
			sprintf(pd->series_labels[3],"T");
			sprintf(pd->series_labels[4],"N");
			for(i = 0; i < 5;i++){
				pd->show_series[i] = 1;
			}
			
			sprintf(pd->plot_title, "Composition of MAPQ >= 20 Reads.");
			sprintf(pd->description,"A HMM was trained on a subset of the sequences. Shown are log2 odds ratios comparing emission probabilities in match states to background nucleotide probabilities. Values above 0 indicate positional enrichment of a particular nucleotide. \"L\" indicates the emission probabilities for a state modelling residiues in the middle of the reads. ");
			//fprintf(stderr,"Got here\n");
			for(j = 2; j <seq_stats->min_len +2;j++){
				
				if(j ==  2 +seq_stats->min_len  / 2){
					sprintf(pd->labels[j-2], "L");
				}else{
					sprintf(pd->labels[j-2], "%d",j-1);
				}
				//	fprintf(stderr,"Got here %d (%d)\n",j ,seq_stats->min_len +2 );
				for(c = 0; c < 5;c++){
					pd->data[c][j-2] =  log2f(scaledprob2prob(hmms[0]->emissions[j][c] ) /( (float) seq_stats->nuc_num[c] / (float)sum));
				}
			}
			//fprintf(stderr,"Got here\n");
			pd->num_points_shown = 20;
			
			pd->num_points = seq_stats->min_len;
			pd->num_series = 4;
			pd->plot_type = LINE_PLOT;
			print_html5_chart(outfile, pd);
			
		}
		
		sum = 0;
		for(i = 0; i < 6;i++){
			sum +=seq_stats->nuc_num[i];
		}
		pd->color_scheme = 0;
		
		if(seq_stats->min_len > 41){
			seq_stats->min_len = 41;
		}
		
		if(hmms[1]){
			
			sprintf(pd->series_labels[0],"A");
			sprintf(pd->series_labels[1],"C");
			sprintf(pd->series_labels[2],"G");
			sprintf(pd->series_labels[3],"T");
			sprintf(pd->series_labels[4],"N");
			for(i = 0; i < 5;i++){
				pd->show_series[i] = 1;
			}
			
			sprintf(pd->plot_title, "Composition of  0 >= MAPQ <  20 Reads. ");
			sprintf(pd->description,"A HMM was trained on a subset of the sequences. Shown are log2 odds ratios comparing emission probabilities in match states to background nucleotide probabilities. Values above 0 indicate positional enrichment of a particular nucleotide. \"L\" indicates the emission probabilities for a state modelling residiues in the middle of the reads. ");
			for(j = 2; j <seq_stats->min_len +2;j++){
				if(j ==  2 +seq_stats->min_len  / 2){
					sprintf(pd->labels[j-2], "L");
				}else{
					sprintf(pd->labels[j-2], "%d",j-1);
				}
				for(c = 0; c < 5;c++){
					pd->data[c][j-2] =  log2f(scaledprob2prob(hmms[1]->emissions[j][c] ) /( (float) seq_stats->nuc_num[c] / (float)sum));
				}
			}
			pd->num_points_shown = 20;
			
			pd->num_points = seq_stats->min_len;
			pd->num_series = 4;
			pd->plot_type = LINE_PLOT;
			print_html5_chart(outfile, pd);
		}
		
		pd->color_scheme = 0;
		if(hmms[2]){
			sprintf(pd->series_labels[0],"A");
			sprintf(pd->series_labels[1],"C");
			sprintf(pd->series_labels[2],"G");
			sprintf(pd->series_labels[3],"T");
			sprintf(pd->series_labels[4],"N");
			for(i = 0; i < 5;i++){
				pd->show_series[i] = 1;
			}
			
			sprintf(pd->plot_title, "Composition of unmapped reads.");
			sprintf(pd->description,"A HMM was trained on a subset of the sequences. Shown are log2 odds ratios comparing emission probabilities in match states to background nucleotide probabilities. Values above 0 indicate positional enrichment of a particular nucleotide. \"L\" indicates the emission probabilities for a state modelling residiues in the middle of the reads.");
			for(j = 2; j <seq_stats->min_len +2;j++){
				if(j ==  2 +seq_stats->min_len  / 2){
					sprintf(pd->labels[j-2], "L");
				}else{
					sprintf(pd->labels[j-2], "%d",j-1);
				}
				for(c = 0; c < 5;c++){
					pd->data[c][j-2] =  log2f(scaledprob2prob(hmms[2]->emissions[j][c] ) /( (float) seq_stats->nuc_num[c] / (float)sum));
				}
			}
			pd->num_points_shown = 20;
			
			pd->num_points = seq_stats->min_len;
			pd->num_series = 4;
			pd->plot_type = LINE_PLOT;
			print_html5_chart(outfile, pd);
			
		}
		
		
		
		
		
		sprintf(pd->series_labels[0],"A");
		sprintf(pd->series_labels[1],"C");
		sprintf(pd->series_labels[2],"G");
		sprintf(pd->series_labels[3],"T");
		sprintf(pd->series_labels[4],"N");
		
		for(i = 0; i < 6;i++){
			if(!seq_stats->alignments[i]){
				pd->show_series[i] =0;
			}
		}
		
		for(i = 0; i < 5;i++){
			switch (i) {
				case 0:
					sprintf(pd->plot_title, "Distribution of Mismatches (MAPQ >= 30):");
					sprintf(pd->description,"Distribution of Mismatches in MAPQ >= 30 reads.");
					break;
				case 1:
					sprintf(pd->plot_title, "Distribution of Mismatches (MAPQ < 30):");
					sprintf(pd->description,"Distribution of Mismatches in MAPQ < 30 reads.");
					break;
					
				case 2:
					sprintf(pd->plot_title, "Distribution of Mismatches (MAPQ < 20):");
					sprintf(pd->description,"Distribution of Mismatches in MAPQ < 20 reads.");
					break;
					
				case 3:
					sprintf(pd->plot_title, "Distribution of Mismatches (MAPQ < 10):");
					sprintf(pd->description,"Distribution of Mismatches in MAPQ < 10 reads.");
					break;
					
				case 4:
					sprintf(pd->plot_title, "Distribution of Mismatches(MAPQ < 3):");
					sprintf(pd->description,"Distribution of Mismatches in MAPQ < 3 reads.");
					break;
					
				default:
					break;
			}
			if(seq_stats->alignments[i]){
				for(j = 0; j <= seq_stats->max_len;j++){
					sprintf(pd->labels[j], "%dnt",j+1);
					for(c = 0; c < 5;c++){
						pd->data[c][j] =  (float)seq_stats->mismatches[i][j][c] / (float)seq_stats->alignments[i] * 100.0f;
					}
				}
			}
			pd->num_points = seq_stats->max_len;
			pd->num_series = 4;
			pd->plot_type = BAR_PLOT;
			print_html5_chart(outfile, pd);
		}
		
		pd->width = 900;
		pd->color_scheme = 0;
		
		//fprintf(stderr,"Errors :\n");
		for(i = 0; i < 5;i++){
			switch (i) {
				case 0:
					sprintf(pd->plot_title, "Number of Errors Per Read (MAPQ >= 30):");
					sprintf(pd->description,"Barplot shows the percentage of reads (y-axis) with 0, 1, 2 ... errors (x axis) for MAPQ >= 30 reads.");
					break;
				case 1:
					sprintf(pd->plot_title, "Number of Errors Per Read(MAPQ < 30):");
					sprintf(pd->description,"Barplot shows the percentage of reads (y-axis) with 0, 1, 2 ... errors (x axis) for MAPQ < 30 reads.");
					break;
					
				case 2:
					sprintf(pd->plot_title, "Number of Errors Per Read(MAPQ < 20):");
					sprintf(pd->description,"Barplot shows the percentage of reads (y-axis) with 0, 1, 2 ... errors (x axis) for MAPQ < 20 reads.");
					break;
					
				case 3:
					sprintf(pd->plot_title, "Number of Errors Per Read(MAPQ < 10):");
					sprintf(pd->description,"Barplot shows the percentage of reads (y-axis) with 0, 1, 2 ... errors (x axis) for MAPQ < 10 reads.");
					break;
					
				case 4:
					sprintf(pd->plot_title, "Number of Errors Per Read(MAPQ < 3):");
					sprintf(pd->description,"Barplot shows the percentage of reads (y-axis) with 0, 1, 2 ... errors (x axis) for MAPQ < 3 reads.");
					break;
					
				default:
					break;
			}
			if(seq_stats->alignments[i]){
				for(j = 0; j <= seq_stats->max_error_per_read;j++){
					sprintf(pd->labels[j], "%d",j);
					pd->data[0][j] = 100.0 * (float)seq_stats->errors[i][j]/ (float)seq_stats->alignments[i];
				}
				
				sprintf(pd->series_labels[0],"Errors");
				pd->num_points = seq_stats->max_error_per_read ;
				pd->num_series = 1;
				pd->width = 300;
				pd->num_points_shown = 10;//seq_stats->max_error_per_read;
				
				
				pd->plot_type = BAR_PLOT;
				print_html5_chart(outfile, pd);
			}
		}
		
		for(i =0 ; i < 3;i++){
			if(hmms[i]){
				free_hmm(hmms[i]);
			}
		}
		MFREE(hmms);
		
		print_html5_footer(outfile);
		
		free_plot_data(pd);
		
		fclose(outfile);
		sprintf(param->buffer,"\n\n");
		param->messages = append_message(param->messages, param->buffer);
		
	}
	
	
	MFREE(hmm_data->length);
	MFREE(hmm_data->score);
	MFREE(hmm_data->string);
	MFREE(hmm_data->weight);
	MFREE(hmm_data);

	///
	
	
	free_seq_stats(seq_stats);
	
	free_read_info(ri, param->num_query);
	free_param(param);
	
	
	return EXIT_SUCCESS;
	
}

char* make_file_stats(char* filename,char* buffer)
{
	struct stat buf;
	int local_ret;
	
	char time_string[200];
	int hour;
	struct tm *ptr;
		
	local_ret= stat ( filename, &buf );
	if ( local_ret == 0 ){
		// %9jd", (intmax_t)statbuf.st_size);
		
		ptr = localtime(&buf.st_mtime);
		hour = ptr->tm_hour;
		if (hour <= 11){
		//	am_or_pm = 'a';
		}else {
			hour -= 12;
		//	am_or_pm = 'p';
		}
		if (hour == 0){
			hour = 12;
		}
		strftime(time_string, 200, "%F %H:%M:%S\t", ptr);

		sprintf(buffer,"size:%9d bytes, created %s",(int)buf.st_size, time_string);
	}else{
		fprintf(stderr,"Failed getting stats for file:%s\n",filename );
	}
	
	
	return buffer;
}


struct seq_stats* reformat_base_qualities(struct seq_stats* seq_stats)
{
	//From wikipedia:
	/*SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
	..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
	...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
	.................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
	LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
	!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
	|                         |    |        |                              |                     |
	33                        59   64       73                            104                   126
	0........................26...31.......40
	-5....0........9.............................40
	0........9.............................40
	3.....9.............................40
	0.2......................26...31.........41
	
	S - Sanger        Phred+33,  raw reads typically (0, 40)
	X - Solexa        Solexa+64, raw reads typically (-5, 40)
	I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
	J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
	with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
	(Note: See discussion above).
	L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)*/
	int start = -1;
	int stop = -1;
	//fprintf(stderr,"Got here\n");
	int i;
	for(i =0;i < 256;i++){
		if(seq_stats->base_qualities[i]){
			if(start == -1){
				start = i;
			}
			stop = i;
		}
	}

	if(start == 33 && stop == 73 ){
#ifdef DEBUG
		fprintf(stderr,"S - Sanger\n");
#endif
	}
	if(start == 49 && stop == 104){
#ifdef DEBUG
		fprintf(stderr,"X - SOLEXA\n");
#endif
	}
	if(start == 64 && stop == 104){
#ifdef DEBUG
		fprintf(stderr,"Illumina 1.3+\n");
#endif
	}
	if(start == 66&& stop == 104 ){
#ifdef DEBUG
		fprintf(stderr,"Illumina 1.5+\n");
#endif
	}
	if(start == 33 && stop == 74){
#ifdef DEBUG
		fprintf(stderr,"Illumina 1.8+\n");
#endif
	}

	seq_stats->base_quality_offset = start;
	return seq_stats;
}




struct hmm* init_samstat_hmm(int average_length, int max_sequence_len)
{
	struct hmm* hmm = NULL;
	
	if(average_length > 41){
		average_length = 41;
	}
	
	hmm = malloc_hmm(average_length + 2, 5,max_sequence_len+2);
	
	init_logsum();
	
	int i,j,c;
	for(i = 0; i < hmm->num_states;i++){
		for(j = 0; j < hmm->num_states;j++){
			hmm->transitions[i][j] = prob2scaledprob(0.0);
			hmm->transitions_e[i][j] = prob2scaledprob(0.5);
		}
	}
	
	hmm->transitions[STARTSTATE ][2]= prob2scaledprob(0.5f);
	//hmm->transitions[STARTSTATE][3] = prob2scaledprob(0.5f);
	c = 2;
	for(i = 0; i < average_length / 2 ;i++){
		hmm->transitions[c][c+1]= prob2scaledprob(0.5f);
		c++;
		//hmm->transitions[i+2][ENDSTATE]= prob2scaledprob(0.5f);
	}
	hmm->transitions[c][c]= prob2scaledprob(0.5f);
	hmm->transitions[c][c+1]= prob2scaledprob(0.5f);
	//hmm->transitions[i][ENDSTATE]= prob2scaledprob(0.5f);
	for(i =c+1 ; i < average_length+2  ;i++){
		hmm->transitions[c][c+1]= prob2scaledprob(0.5f);
		c++;
		//hmm->transitions[i+2][ENDSTATE]= prob2scaledprob(0.5f);
	}
	
	hmm->transitions[c][ENDSTATE]= prob2scaledprob(0.5f);

	
	
	//norm;
	float sum = prob2scaledprob(0.0);
	for(i = 0; i < hmm->num_states;i++){
		if(i != 1){
			sum = prob2scaledprob(0.0);
			for(j = 0; j < hmm->num_states;j++){
				sum =  logsum(sum ,hmm->transitions[i][j] ) ;
				if(hmm->transitions[i][j] == -INFINITY){
					hmm->transitions_e[i][j] = prob2scaledprob(0.0f);
				}
			}
			
			for(j = 0; j < hmm->num_states;j++){
				hmm->transitions[i][j] = hmm->transitions[i][j] -sum;
			}
		}
	}
	
	
	
	for(i =1 ; i < hmm->num_states;i++){
		MMALLOC(hmm->emissions[i], sizeof(float) * hmm->alphabet_len);
		MMALLOC(hmm->emissions_e[i], sizeof(float) * hmm->alphabet_len);
		for(j = 0;j < hmm->alphabet_len;j++){
			hmm->emissions[i][j] = prob2scaledprob(1.0f / (float) hmm->alphabet_len);
			hmm->emissions_e[i][j] = prob2scaledprob(0.5);
		}
	}
	
	for(i = 0; i < hmm->num_states;i++){
		c = 0;
		for(j = 0; j < hmm->num_states;j++){
			if(hmm->transitions[i][j]  != -INFINITY){
				hmm->tindex[i][c+1] = j;
				c++;
			}
		}
		hmm->tindex[i][0] = c+1;
	}
	
	//print_hmm_parameters(hmm);
	//exit(-1);
	return hmm;

}

struct seq_stats* init_seq_stats(void)
{
	struct seq_stats* seq_stats = NULL;
	int i,j,c;
	
	MMALLOC(seq_stats, sizeof(struct seq_stats));
	
	seq_stats->alignments = NULL;
	seq_stats->aln_quality = NULL;
	seq_stats->deletions = NULL;
	seq_stats->errors = NULL;
	seq_stats->insertions =NULL;
	seq_stats->mismatches = NULL;
	seq_stats->nuc_composition = NULL;
	seq_stats->nuc_num = NULL;
	seq_stats->percent_identity = NULL;
	seq_stats->seq_len = NULL;
	seq_stats->seq_quality = NULL;
	seq_stats->seq_quality_count = NULL;
	
	seq_stats->base_qualities = NULL;
	seq_stats->total_reads = 0;
	
	MMALLOC(seq_stats->base_qualities, sizeof(int)* 256);
	MMALLOC(seq_stats->alignments, sizeof(int)* 6);
	
	MMALLOC(seq_stats->seq_len,sizeof(int*)* 6);
	MMALLOC(seq_stats->nuc_composition,sizeof(int**)* 6);
	MMALLOC(seq_stats->seq_quality,sizeof(int*)* 6);
	MMALLOC(seq_stats->seq_quality_count,sizeof(int*)* 6);
	MMALLOC(seq_stats->aln_quality,sizeof(int)*6);
	MMALLOC(seq_stats->nuc_num,sizeof(int) * 6);
	//seq_stats->overall_kmers= malloc(sizeof(float) * KMERALLOC);
	MMALLOC(seq_stats->errors,sizeof(float*)* 6 );
	
	MMALLOC(seq_stats->percent_identity,sizeof(double)* 6 );
	
	MMALLOC(seq_stats->mismatches,sizeof(int**)* 6);
	MMALLOC(seq_stats->insertions, sizeof(int**) * 6);
	MMALLOC(seq_stats->deletions, sizeof(int*) * 6);
	
	for(i= 0; i < 256;i++){
		seq_stats->base_qualities[i] = 0;
	}
	
	for(c = 0; c < 6;c++){
		seq_stats->mismatches[c] = NULL;
		seq_stats->insertions[c] = NULL;
		seq_stats->deletions[c] = NULL;
		seq_stats->errors[c] = NULL;
		
		seq_stats->nuc_num[c] = 0;
		seq_stats->alignments[c] = 0;
		
		
		seq_stats->seq_len[c] = NULL;
		seq_stats->nuc_composition[c] = NULL;
		seq_stats->seq_quality[c] = NULL;
		seq_stats->seq_quality_count[c] = NULL;
	
		
		
		
		MMALLOC(seq_stats->mismatches[c],sizeof(int*)* MAX_SEQ_LEN);
		MMALLOC(seq_stats->insertions[c],sizeof(int*) * MAX_SEQ_LEN);
		MMALLOC(seq_stats->deletions[c],sizeof(int) * MAX_SEQ_LEN);
		MMALLOC(seq_stats->errors[c],sizeof(float)* MAXERROR );
		
		MMALLOC(seq_stats->seq_len[c],sizeof(int)* MAX_SEQ_LEN);
		MMALLOC(seq_stats->nuc_composition[c],sizeof(int*)* MAX_SEQ_LEN);
		MMALLOC(seq_stats->seq_quality[c],sizeof(int)* MAX_SEQ_LEN);
		MMALLOC(seq_stats->seq_quality_count[c],sizeof(int)* MAX_SEQ_LEN);

		
		seq_stats->percent_identity[c] = 0.0f;
		        
		for(i = 0; i < MAXERROR;i++){
			seq_stats->errors[c][i] = 0;
			
		}
		
		for(i= 0; i < MAX_SEQ_LEN;i++){
			seq_stats->deletions[c][i] = 0;
			seq_stats->mismatches[c][i] = NULL;
			seq_stats->insertions[c][i] = NULL;
			seq_stats->seq_len[c][i] = 0;
			seq_stats->nuc_composition[c][i] = NULL;
			
			seq_stats->seq_quality[c][i] = 0;
			seq_stats->seq_quality_count[c][i] = 0;
			
			MMALLOC(seq_stats->mismatches[c][i],sizeof(int)*5);
			MMALLOC(seq_stats->insertions[c][i],sizeof(int)*5);
			MMALLOC(seq_stats->nuc_composition[c][i],sizeof(int*)* 5);
			for(j = 0; j < 5;j++){
				seq_stats->mismatches[c][i][j] = 0;
				seq_stats->insertions[c][i][j] = 0;
				seq_stats->nuc_composition[c][i][j] = 0;
			}
		}
	}
	
	seq_stats->sam = 0;
	seq_stats->md = 0;
	seq_stats->base_quality_offset = 33;
	seq_stats->min_len = 1000000000;
	seq_stats->max_len = -1000000000;
	seq_stats->max_error_per_read =-1000000000;
	return seq_stats;
}

struct seq_stats* clear_seq_stats(struct seq_stats* seq_stats)
{
	//struct seq_stats* seq_stats = NULL;
	int i,j,c;
	
	seq_stats->total_reads = 0;
	
	for(i= 0; i < 256;i++){
		seq_stats->base_qualities[i] = 0;
	}
	
	for(c = 0; c < 6;c++){
		
		seq_stats->nuc_num[c] = 0;
		seq_stats->alignments[c] = 0;
		seq_stats->percent_identity[c] = 0.0f;
		
		for(i = 0; i < MAXERROR;i++){
			seq_stats->errors[c][i] = 0;
			
		}
		
		for(i= 0; i < MAX_SEQ_LEN;i++){
			seq_stats->deletions[c][i] = 0;
			seq_stats->seq_len[c][i] = 0;
			
			seq_stats->seq_quality[c][i] = 0;
			seq_stats->seq_quality_count[c][i] = 0;
			
			for(j = 0; j < 5;j++){
				seq_stats->mismatches[c][i][j] = 0;
				seq_stats->insertions[c][i][j] = 0;
				seq_stats->nuc_composition[c][i][j] = 0;
			}
		}
	}
	
	seq_stats->sam = 0;
	seq_stats->md = 0;
	seq_stats->base_quality_offset = 33;
	seq_stats->min_len = 1000000000;
	seq_stats->max_len = -1000000000;
	seq_stats->max_error_per_read =-1000000000;
	return seq_stats;
}

void free_seq_stats(struct seq_stats* seq_stats)
{
	int i,j;
	
	for(j = 0; j < 6;j++){
		for(i= 0; i < MAX_SEQ_LEN;i++){
			free(seq_stats->mismatches[j][i]);// = malloc(sizeof(int)*5);
			free(seq_stats->insertions[j][i]);// = malloc(sizeof(int)*5);
		}
		free(seq_stats->mismatches[j]);// = malloc(sizeof(int*)* MAX_SEQ_LEN);
		free(seq_stats->insertions[j]);// = malloc(sizeof(int*) * MAX_SEQ_LEN);
		free(seq_stats->deletions[j]);// = malloc(sizeof(int) * MAX_SEQ_LEN);
	}
	free(seq_stats->mismatches);// = malloc(sizeof(int*)* MAX_SEQ_LEN);
	free(seq_stats->insertions);// = malloc(sizeof(int*) * MAX_SEQ_LEN);
	free(seq_stats->deletions);// = malloc(sizeof(int) * MAX_SEQ_LEN);
	
	for(i = 0; i < 6;i++){
		
		for(j = 0; j < MAX_SEQ_LEN;j++){
			free(seq_stats->nuc_composition[i][j]);// = malloc(sizeof(int*)* 5);
			//free(seq_stats->seq_quality[i][j]);// = malloc(sizeof(int)* 6);
			
		}
		//free(seq_stats->seq_quality[i]);
		free(seq_stats->seq_len[i]);// = malloc(sizeof(int)* MAX_SEQ_LEN);
		free(seq_stats->nuc_composition[i]);// = malloc(sizeof(int*)* MAX_SEQ_LEN);
		free(seq_stats->seq_quality[i]);// = malloc(sizeof(int*)* MAX_SEQ_LEN);
		free(seq_stats->seq_quality_count[i]);
	}
	for(i = 0; i < 6;i++){
		free(seq_stats->errors[i]);
		//free(seq_stats->percent_identity[i]);
		
	}
	free(seq_stats->errors);
	free(seq_stats->percent_identity);
	free(seq_stats->base_qualities);
	free(seq_stats->alignments);
	free(seq_stats->nuc_num);
	free(seq_stats->seq_len);// = malloc(sizeof(int*)* 6);
	free(seq_stats->nuc_composition);// = malloc(sizeof(int**)* 6);
	free(seq_stats->seq_quality);// = malloc(sizeof(int**)* 6);
	free(seq_stats->seq_quality_count);
	free(seq_stats->aln_quality);// = malloc(sizeof(int)*6);
	free(seq_stats);// = malloc(sizeof(struct seq_stats));
	
}

 void print_stats(struct seq_stats* seq_stats)
{
	int i,j,c;
	
	fprintf(stderr,"Nucleotides:\n");
	for(c = 0; c < 6;c++){
		fprintf(stderr,"%d	%d\n",c,seq_stats->nuc_num[c]);
	}
	
	fprintf(stderr,"Alignments :\n");
	for(c = 0; c < 6;c++){
		fprintf(stderr,"%d	%d\n",c,seq_stats->alignments[c]);
	}
	
	fprintf(stderr,"Percentage Identitiy :\n");
	for(c = 0; c < 6;c++){
		fprintf(stderr,"%d	%f\n",c,seq_stats->percent_identity[c]);
	}
	
	fprintf(stderr,"Errors :\n");
	for(c = 0; c < 6;c++){
		if(seq_stats->alignments[c]){
		fprintf(stderr,"Class:%d\n",c);
		for(i = 0; i < MAXERROR;i++){
			fprintf(stderr," %f",seq_stats->errors[c][i]);
			
		}
		fprintf(stderr,"\n");
		}
	}
	
	fprintf(stderr,"Quality:  %d\n", seq_stats->base_quality_offset );
	for(c = 0; c < 6;c++){
		if(seq_stats->alignments[c]){
		fprintf(stderr,"Class:%d\n",c);
		for(i= 0; i <= seq_stats->max_len;i++){
			
			fprintf(stderr," %d",seq_stats->seq_quality[c][i] );
			
		}
		fprintf(stderr,"\n");
		}
	}
	
	fprintf(stderr,"Length :\n");
	for(c = 0; c < 6;c++){
		if(seq_stats->alignments[c]){
		fprintf(stderr,"Class:%d\n",c);
		for(i= 0; i <= seq_stats->max_len;i++){
			
			fprintf(stderr," %d",seq_stats->seq_len[c][i] );
			
		}
		fprintf(stderr,"\n");
		}
	}

	
	fprintf(stderr,"Deletions :\n");
	for(c = 0; c < 6;c++){
		if(seq_stats->alignments[c]){
		fprintf(stderr,"Class:%d\n",c);
		for(i= 0; i <= seq_stats->max_len;i++){
		
			fprintf(stderr," %d",seq_stats->deletions[c][i] );
			
		}
		fprintf(stderr,"\n");
		}
	}
	
	
	
	
	
	fprintf(stderr,"Mismatches / Insertions / composition  :\n");
	
	for(c = 0; c < 6;c++){
		if(seq_stats->alignments[c]){
		fprintf(stderr,"Class:%d\n",c);
		
		for(i= 0; i <= seq_stats->max_len;i++){
			fprintf(stderr,"Pos:%d ",i);
			for(j = 0; j < 5;j++){
				fprintf(stderr," %d",seq_stats->mismatches[c][i][j]);
			}
			for(j = 0; j < 5;j++){
				fprintf(stderr," %d",seq_stats->insertions[c][i][j]);
			}
			for(j = 0; j < 5;j++){
				
				fprintf(stderr," %d",seq_stats->nuc_composition[c][i][j]);
				
			}
			fprintf(stderr,"\n");
		
		}
		}
		
	}
}




int parse_cigar_md(struct read_info* ri,struct seq_stats* seq_stats,int qual_key)
{
	int* read = malloc(sizeof(int)* MAX_SEQ_LEN);
	int* genome = malloc(sizeof(int) * MAX_SEQ_LEN);
	int reverse_int[5]  ={3,2,1,0,4};
	char tmp_num[8];
	int l,i,j,c,rp,gp,sp,exit_loop,aln_len,add;
	
	for(i = 0; i < MAX_SEQ_LEN;i++){
		genome[i] = 0;
		read[i] = 0;
	}
	
	l = (int) strlen((char*)ri->cigar);
	exit_loop = 0;
	i =0;
	rp = 0;
	sp = 0;
	while(!exit_loop){
		c = 0;
		if(isdigit((int)ri->cigar[i])){
			j = 0;
			while (isdigit(ri->cigar[i])) {
				tmp_num[j] = ri->cigar[i];
				j++;
				i++;
				if(i == l){
					exit_loop =1;
					break;
				}
			}
			tmp_num[j] = 0;
			
			c = atoi(tmp_num);
		}
		if(isalpha((int)ri->cigar[i])){
			switch (ri->cigar[i]) {
				// Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
					
				case 'M':
				case 'S':
				case '=':
				case 'X':
					for(j = 0; j < c;j++){
						read[rp] = ri->seq[sp];
						rp++;
						sp++;
					}
					//			fprintf(stderr,"M:%d\n",c);
					break;
				case 'I':
					for(j = 0; j < c;j++){
						read[rp] = ri->seq[sp];
						genome[rp] = -1;
						rp++;
						sp++;
					}
					
					//			fprintf(stderr,"I:%d\n",c);
					break;
				case 'D':
					for(j = 0; j < c;j++){
						read[rp] = -1;
						rp++;
						
					}
					
					//			fprintf(stderr,"D:%d\n",c);
					break;
				default:
					break;
			}
			c = 0;
		}
		i++;
		if(i == l){
			exit_loop =1;
			break;
		}
		
	}
	aln_len = rp;
	
	i =0;
	rp = 0;
	
	while(read[rp] == -1){
		rp++;
	}
	gp = 0;
	exit_loop = 0;
	add  = 0;
	l = (int) strlen((char*)ri->md);
	
	//int gg;
	
	//#1C20^ATT0C3A0
	while(!exit_loop){
		if(isdigit((int)ri->md[i])){
			j = 0;
			while (isdigit(ri->md[i])) {
				tmp_num[j] = ri->md[i];
				j++;
				i++;
				if(i == l){
					exit_loop = 1;
					break;
				}
			}
			tmp_num[j] = 0;
			
			c = atoi(tmp_num);
			
			//fprintf(stderr,"MD:%d\n",c);
			for(j = 0; j < c;j++){
				while(genome[gp] == -1){
					gp++;
					rp++;
				}
				while(read[rp] == -1){
					rp++;
				}
				genome[gp] = read[rp];
				
				gp++;
				rp++;
				//fprintf(stderr,"%d	%d	%d \n",aln_len,rp,i);
				while(read[rp] == -1){
					rp++;
				}
			}
			add = 0;
		}else if(isalpha((int)ri->md[i])){
			//fprintf(stderr,"MD:%c\n",ri->md[i]);
			while(genome[gp] == -1){
				gp++;
			}
			genome[gp] = nuc_code[(int)ri->md[i]];
			gp++;
			i++;
			if(!add){
				rp++;
			}
			
		}else{
			add = 1;
			i++;
		}
		if(i == l){
			exit_loop = 1;
			break;
		}
	}
	
	if(ri->strand == 0){
		gp = 0;
		for(i =0;i < aln_len;i++){
			
			if(read[i] != -1 && genome[i] != -1){
				if(read[i] != genome[i]){
					seq_stats->mismatches[qual_key][gp][read[i]] += 1;
				}
				gp++;
				//			fprintf(stderr,"Mismatch %d\n",i);
			}else if(read[i] == -1 && genome[i] != -1){
				seq_stats->deletions[qual_key][gp] += 1;
				
				//			fprintf(stderr,"Deletion %d\n",i);
			}else if(read[i] != -1 && genome[i] == -1){
				seq_stats->insertions[qual_key][gp][read[i]] += 1;
				gp++;
				//			fprintf(stderr,"Insertion %d\n",i);
			}
		}
	}else{
		gp = ri->len-1;
		for(i = 0;i < aln_len;i++){
			
			if(read[i] != -1 && genome[i] != -1){
				if(read[i] != genome[i]){
					seq_stats->mismatches[qual_key][gp][reverse_int[read[i]]] += 1;
					//		fprintf(stderr,"Mismatch %d->%d\n",i,gp);
				}
				gp--;
				
			}else if(read[i] == -1 && genome[i] != -1){
				seq_stats->deletions[qual_key][gp] += 1;
				
				//	fprintf(stderr,"Deletion %d\n",i);
			}else if(read[i] != -1 && genome[i] == -1){
				seq_stats->insertions[qual_key][gp][reverse_int[read[i]]] += 1;
				gp--;
				//	fprintf(stderr,"Insertion %d\n",i);
			}
		}
	}
	
	free(read);
	free(genome);
	return aln_len;
}











