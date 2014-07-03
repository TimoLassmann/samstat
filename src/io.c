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

/*! \file io.c
 \brief functions for reading sequences.
 
 Initializes nucleotide alphabet needed to parse input. Calls parameter parser. Calls functions to process the data. \author Timo Lassmann \bug No known bugs.
 */

#include <ctype.h>
#include "interface.h"
#include "nuc_code.h"
#include "misc.h"

#include "io.h"

#ifndef MMALLOC
#include "malloc_macro.h"
#endif



/** \fn int qsort_ri_prob_compare(const void *a, const void *b)
 \brief Compares reads based their probability.
 Used to sort arrays of string using qsort.
 \param a void pointer to first @ref read_info.
 \param b void pointer to second @ref read_info.
 */
int qsort_ri_mapq_compare(const void *a, const void *b)
{
	
	//struct mys **a = (struct mys **)i1;
	//struct mys **b = (struct mys **)i2;
	//return (*b)->id - (*a)->id;
	
	const struct read_info **elem1 = (const struct read_info**) a;
	
	const struct read_info **elem2 = (const struct read_info**) b;
	
	if ( (*elem1)->mapq > (*elem2)->mapq)
		return -1;
	
	else if ((*elem1)->mapq < (*elem2)->mapq)
		return 1;
	
	else
		return 0;
}



/** \fn FILE* io_handler(FILE* file, int file_num,struct parameters* param)
 \brief Opens stream to files. 
 
 Recognizes file types by prefix. 
 
 Used to sort arrays of string using qsort.
 \param file empty file pointer.
\param file_num index of input file.
 \param param @ref parameters.

 */

FILE* io_handler(FILE* file, int file_num,struct parameters* param)
{
	char command[1000];
	char  tmp[1000];
	int i = 0; 
	int gzcat = -1;
	if(access("/usr/bin/gzcat", X_OK) == 0){
		gzcat = 1;
	}else if(access("/bin/gzcat", X_OK) == 0){
		gzcat = 1;
	}else if(access("/usr/bin/zcat", X_OK) == 0){
		gzcat = 0;
	}else if(access("/bin/zcat", X_OK) == 0){
		gzcat = 0;
	}

	param->gzipped = 0;
	param->bzipped = 0;
	param->sam = 0;
	param->fasta = 0;
	
	if(!file_exists(param->infile[file_num])){
		sprintf(param->buffer,"Error: Cannot find input file: %s\n",param->infile[file_num] );
		param->messages = append_message(param->messages, param->buffer);
		free_param(param);
		exit(EXIT_FAILURE);
	}
	
	if(!strcmp(".sam", param->infile[file_num] + (strlen(param->infile[file_num] ) - 4))){
		param->sam = 1;
	}else if (!strcmp(".bam", param->infile[file_num] + (strlen(param->infile[file_num] ) - 4))){
		param->sam = 2;
	}else if (!strcmp(".fa", param->infile[file_num] + (strlen(param->infile[file_num] ) - 3))){
		param->sam = 0;
		param->fasta = 1;
	}else if (!strcmp(".fq", param->infile[file_num] + (strlen(param->infile[file_num] ) - 3))){
		param->sam = 0;
	}else if (!strcmp(".fastq", param->infile[file_num] + (strlen(param->infile[file_num] ) - 6))){
		param->sam = 0;
	}else if (!strcmp(".fastaq", param->infile[file_num] + (strlen(param->infile[file_num] ) - 7))){
		param->sam = 0;
	}else if (!strcmp(".fasta", param->infile[file_num] + (strlen(param->infile[file_num] ) - 6))){
		param->sam = 0;
		param->fasta = 1;
	}else if(!strcmp(".sam.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 7))){
		param->sam = 1;
		param->gzipped  = 1;
	}else if (!strcmp(".bam.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 7))){
		param->sam = 2;
		param->gzipped  = 1;
	}else if (!strcmp(".fa.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 6))){
		param->sam = 0;
		param->fasta = 1;
		param->gzipped  = 1;
	}else if (!strcmp(".fq.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 6))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fastq.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 9))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fastaq.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 10))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fasta.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 9))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fastq.bz2", param->infile[file_num] + (strlen(param->infile[file_num] ) - 10))){
		param->sam = 0;
		param->bzipped  = 1;
	}else if (!strcmp(".fq.bz2", param->infile[file_num] + (strlen(param->infile[file_num] ) - 7))){
		param->sam = 0;
		param->bzipped  = 1;
	}else{
		param->sam = -1;
	}
	
	
	if(param->gzipped && gzcat == -1){
		sprintf(param->buffer,"Cannot find gzcat / zcat on your system. Try gzcat <infile> | samstat -f sam/bam/fa/fq\n");
		param->messages = append_message(param->messages, param->buffer);
		free_param(param);
		exit(EXIT_FAILURE);
	}
	
	if(file_num == -1){
		if(param->sam == 2){
			command[0] = 0;
			if(!param->filter){
				strcat ( command, "samtools view -F 768 "); 
			}else{
				strcat ( command, "samtools view -F "); 
				i = sprintf (tmp, "%s ",param->filter);
				strcat ( command, tmp);
			}
			i = sprintf (tmp, "%s ","-");
			strcat ( command, tmp);
			if (!(file = popen(command, "r"))) {
				fprintf(stderr,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				exit(-1);
			}
		}else if(param->sam == 1){
			command[0] = 0;
			if(!param->filter){
				strcat ( command, "samtools view -SF 768 "); 
			}else{
				strcat ( command, "samtools view -SF "); 
				i = sprintf (tmp, "%s ",param->filter);
				strcat ( command, tmp);
			}
			i = sprintf (tmp, "%s ", "-");
			strcat ( command, tmp);
			if (!(file = popen(command, "r"))) {
				fprintf(stderr,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				exit(-1);
			}
		}else{
			file = stdin;
		}
	}else{
		if(param->sam == 2){
			command[0] = 0;
			
			if(param->bzipped){
				strcat ( command, "bzcat ");
				if(!param->filter){
					i = sprintf (tmp, "%s | samtools view -F 768 - ", param->infile[file_num]);
					strcat ( command, tmp);
				}else{
					i = sprintf (tmp, "%s | samtools view -F  ", param->infile[file_num]);
					strcat ( command, tmp);
					i = sprintf (tmp, "%s - ",param->filter);
					strcat ( command, tmp);
				}
				
			}else if(param->gzipped){
				if(gzcat == 1){
					strcat ( command, "gzcat "); 
				}else{
					strcat ( command, "zcat "); 
				}
				if(!param->filter){
					i = sprintf (tmp, "%s | samtools view -F 768 - ", param->infile[file_num]);
					strcat ( command, tmp);
				}else{
					i = sprintf (tmp, "%s | samtools view -F  ", param->infile[file_num]);
					strcat ( command, tmp);
					i = sprintf (tmp, "%s - ",param->filter);
					strcat ( command, tmp);
				}
			}else{
				if(!param->filter){
					strcat ( command, "samtools view -F 768 "); 
				}else{
					strcat ( command, "samtools view -F "); 
					i = sprintf (tmp, "%s ",param->filter);
					strcat ( command, tmp);
				}
				i = sprintf (tmp, "%s ", param->infile[file_num]);
				strcat ( command, tmp);
			}
			if (!(file = popen(command, "r"))) {
				sprintf(param->buffer,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				param->messages = append_message(param->messages, param->buffer);
				free_param(param);
				exit(EXIT_FAILURE);
			}
		}else if(param->sam == 1){
			command[0] = 0;
			if(param->gzipped){
				if(gzcat == 1){
					strcat ( command, "gzcat "); 
				}else{
					strcat ( command, "zcat "); 
				}
				if(!param->filter){
					i = sprintf (tmp, "%s | samtools view -SF 768 - ", param->infile[file_num]);
					strcat ( command, tmp);
				}else{
					i = sprintf (tmp, "%s | samtools view -SF  ", param->infile[file_num]);
					strcat ( command, tmp);
					i = sprintf (tmp, "%s - ",param->filter);
					strcat ( command, tmp);
				}
			}else{
				if(!param->filter){
					strcat ( command, "samtools view -SF 768 "); 
				}else{
					strcat ( command, "samtools view -SF "); 
					i = sprintf (tmp, "%s ",param->filter);
					strcat ( command, tmp);
				}
				i = sprintf (tmp, "%s ", param->infile[file_num]);
				strcat ( command, tmp);
			}
			if (!(file = popen(command, "r"))) {
				sprintf(param->buffer,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				param->messages = append_message(param->messages, param->buffer);
				free_param(param);
				exit(EXIT_FAILURE);
			}
		}else{
			command[0] = 0;
			if(param->bzipped){
				strcat ( command, "bzcat ");
				
			}else if(param->gzipped){
				if(gzcat == 1){
					strcat ( command, "gzcat ");
				}else{
					strcat ( command, "zcat ");
				}
			}else{
				strcat ( command, "cat ");
			}
			i = sprintf (tmp, "%s ", param->infile[file_num]);
			strcat ( command, tmp);
			//fprintf(stderr,"%s\n",command);
			if (!(file = popen(command, "r"))) {
				sprintf(param->buffer,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				param->messages = append_message(param->messages, param->buffer);
				free_param(param);
				exit(EXIT_FAILURE);
			}
		}
	}
	return file;
}




int read_sam_chunk(struct read_info** ri,struct parameters* param,FILE* file)
{
	char line[MAX_LINE];
	int column = 0; 
	int i,j,g,tmp;
	unsigned int pos;
	int read = 0;
	int c = 0;
	int hit = 0;
	
	ri = clear_read_info(ri, param->num_query);
	
	while(fgets(line, MAX_LINE, file)){
		if(line[0] != '@'){
			column = 1; //<QNAME>
			tmp = 0;
			hit = 0;
			pos = 0xFFFFFFFFu;
			for(j = 0;j < MAX_LINE;j++){
				tmp++;
				if(isspace((int)line[j])){
					break;
				}
			}
			
			MMALLOC(ri[c]->name,sizeof(unsigned char)* tmp);
			for(j = 0;j < MAX_LINE;j++){
				
				if(isspace((int)line[j])){
					ri[c]->name[j] = 0;
					break;
				}
				ri[c]->name[j] = line[j];
			}
			
			for(i = 0; i < MAX_LINE;i++){
				if(line[i] == '\n'){
					break;
				}
				if(isspace((int)line[i])){
					column++;
					switch(column){
						case 2: // <FLAG>
							tmp = atoi(line+i+1);
							ri[i]->strand = (tmp & 0x10);

							//WARNING - read should be reverse complemented if mapped to negative strand before tagdusting...
							
							/*tmp = atoi(line+i+1);
							ri[c]->strand[hit] = (tmp & 0x10);
							if(tmp == 4){
								ri[c]->hits[hit] = 0;
							}else{
								ri[c]->hits[hit] = 1;
							}
							hit++;*/
							
							break;
						case 3: // <RNAME> 
							
							break;
						case 4: // <POS>
							
							break;
						case 5: //  <MAPQ>
							
							ri[c]->mapq =  atof(line +i +1); 
							
							break;
						case 6: //  <CIGAR>
							tmp = 0;
							for(j = i+1;j < MAX_LINE;j++){
								tmp++;
								if(isspace((int)line[j])){
									break;
								}
							}
							
							ri[c]->cigar = malloc(sizeof(unsigned char)* tmp);
							g = 0;
							for(j = i+1;j < MAX_LINE;j++){
								if(isspace((int)line[j])){
									ri[c]->cigar[g] = 0;
									break;
								}
								ri[c]->cigar[g] = line[j];
								g++;
							}
							break;
						case 7: //  <MRNM>
							break;
						case 8: //  <MPOS>
							break;
						case 9: //  <ISIZE>
							break;
						case 10: // <SEQ>
							
							tmp = 0;
							for(j = i+1;j < MAX_LINE;j++){
								tmp++;
								if(isspace((int)line[j])){
									break;
								}
							}
							
							MMALLOC(ri[c]->seq,sizeof(unsigned char)* tmp);
							MMALLOC(ri[c]->labels,sizeof(unsigned char)* tmp);
							
							g = 0;
							for(j = i+1;j < MAX_LINE;j++){
								
								if(isspace((int)line[j])){
									ri[c]->seq[g] = 0;
									ri[c]->labels[g] = 0;
									break;
								}
								ri[c]->seq[g] = nuc_code[(int)line[j]];
								ri[c]->labels[g] = 0;

								g++;
							}
							
							ri[c]->len = g;
							break;
						case 11: // <QUAL>
							tmp = 0;
							for(j = i+1;j < MAX_LINE;j++){
								tmp++;
								if(isspace((int)line[j])){
									break;
								}
							}
							g= 0;
							MMALLOC(ri[c]->qual,sizeof(unsigned char)* tmp);
							for(j = i+1;j < MAX_LINE;j++){
								
								if(isspace((int)line[j])){
									ri[c]->qual[g] = 0;
									break;
								}
								ri[c]->qual[g] = line[j];
								g++;
							}
							break;
						default: 
							
									
							i = MAX_LINE;
							break;
					}				}

			}
			tmp = byg_end("NM:i:", line  );
			if(tmp){
				ri[c]->errors = atoi(line+tmp);
				//if(ri[c]->errors > 20){
				///fprintf(stderr,"%s\n,%c,%c,%c,%d\n",line, *(line +tmp), *(line +tmp+1),*(line +tmp+2), ri[c]->errors);
				//}
				
			}else{
				ri[c]->errors = -1;
			}
			tmp = byg_end("MD:Z:", line  );
			if(tmp){
				g = 0;
				for(j = tmp ;j < MAX_LINE;j++){
					g++;
					if(isspace((int)line[j])){
						break;
					}
					
				}
				ri[c]->md = malloc(sizeof(unsigned char)* g);
				g = 0;
				for(j = tmp ;j < MAX_LINE;j++){
					
					if(isspace((int)line[j])){
						ri[c]->md[g] = 0;
						break;
					}
					ri[c]->md[g] = line[j];
					g++;
				}
			}
						
			
			
			//ri[c]->hits[hit] = 0xFFFFFFFFu;
			
			c++;
			read++;
			if(c == param->num_query){
				return c;
			}
		}
	}
	return c;
}

int read_fasta_fastq(struct read_info** ri,struct parameters* param,FILE *file) 
{
	int park_pos = -1;
	char line[MAX_LINE];
	int i;//,j;
	int seq_p = 0;
	int set = 0;
	int len = 0;
	int size = 0;
	
	ri = clear_read_info(ri, param->num_query);
	while(fgets(line, MAX_LINE, file)){
		if((line[0] == '@' && !set)|| (line[0] == '>' && !set)){
			//set sequence length of previous read
			
			//check if there is still space....
			//if(param->num_query == size){
			//	fseek (file , -  strlen(line) , SEEK_CUR);
			//	return size;
			//}
			park_pos++;
			len = 0;
			seq_p = 1;
			for(i = 1;i < MAX_LINE;i++){
				len++;
				if(iscntrl((int)line[i])){
					break;
				}
				
			}
			
			//ri[park_pos]->hits[0] = 0;
			//ri[park_pos]->strand[0] = 0;
			MMALLOC(ri[park_pos]->name,sizeof(unsigned char)* (len+1));
			for(i = 1;i < MAX_LINE;i++){
				
				if(iscntrl((int)line[i])){
					ri[park_pos]->name[i-1] = 0;
					break;
				}
				if(isspace((int)line[i])){
					ri[park_pos]->name[i-1] = ';';
				}
				
				ri[park_pos]->name[i-1] = line[i];
			}
			//fprintf(stderr,"LEN:%d	%s\n",len,ri[park_pos]->name);
			
			set = 1;
			size++;
			//get ready to read quality if present  
		}else if(line[0] == '+' && !set){
			seq_p = 0;
			set = 1;
			//reading sequence or quality  
		}else{	
			if(set){
				if(seq_p){
					len = 0;
					for(i = 0;i < MAX_LINE;i++){
						len++;
						if(iscntrl((int)line[i])){
							break;
						}
					}
					//fprintf(stderr,"SEQ LEN:%d	%s\n",len,line);
					MMALLOC(ri[park_pos]->seq,sizeof(unsigned char)* (len+1));
					
					MMALLOC(ri[park_pos]->labels, sizeof(unsigned char)* (len+1));
					
					for(i = 0;i < MAX_LINE;i++){
						if(iscntrl((int)line[i])){
							ri[park_pos]->seq[i] = 0;
							ri[park_pos]->labels[i] = 0;
							break;
						}
						ri[park_pos]->seq[i] = nuc_code[(int)line[i]];
						ri[park_pos]->labels[i] = 0;
					}
					ri[park_pos]->len = len-1;
				}else{
					len = 0;
					for(i = 0;i < MAX_LINE;i++){
						len++;
						if(iscntrl((int)line[i])){
							break;
						}
						
					}
					
					if(len-1 != ri[park_pos]->len ){
						sprintf(param->buffer,"ERROR: Length of sequence and base qualities differ!.\n");
						param->messages = append_message(param->messages, param->buffer);
						free_param(param);
						exit(EXIT_FAILURE);
					}
					
					//fprintf(stderr,"QUAL LEN:%d\n",len);
					MMALLOC(ri[park_pos]->qual,sizeof(unsigned char)* (len+1));
					for(i = 0;i < MAX_LINE;i++){
						if(iscntrl((int)line[i])){
							ri[park_pos]->qual[i] = 0;
							break;
						}
						ri[park_pos]->qual[i] = line[i];
					}
				}
			}
			set = 0;
		}
		if(param->num_query == size ){//here I know I am in the last entry AND filled the quality...
			if(!param->fasta && ri[park_pos]->qual){
				return size;
			}
			if(param->fasta && ri[park_pos]->seq){
			   
				return size;
			}
		}
	}
	return size;
}


struct read_info** malloc_read_info(struct read_info** ri, int numseq)
{
	int i;
	MMALLOC(ri, sizeof(struct read_info*) * numseq);
	
	for(i = 0; i < numseq;i++){
		ri[i] = 0;
		MMALLOC(ri[i], sizeof(struct read_info));
		
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->labels = 0;
		ri[i]->len = 0;
		ri[i]->mapq = -1.0;
		ri[i]->cigar = 0;
		ri[i]->md = 0;

		ri[i]->errors = 0;
		ri[i]->strand = 0;
	}
	return ri;
}

struct read_info** clear_read_info(struct read_info** ri, int numseq)
{
	int i;
	
	for(i = 0; i < numseq;i++){
		if(ri[i]->md){
			MFREE(ri[i]->md);
		}

		if(ri[i]->cigar){
			MFREE(ri[i]->cigar);
		}
		if(ri[i]->seq){
			MFREE(ri[i]->seq);
		}
		
		
		
		if(ri[i]->name){
			MFREE(ri[i]->name);
		}
		if(ri[i]->qual){
			MFREE(ri[i]->qual);
		}
		if(ri[i]->labels){
			MFREE(ri[i]->labels);
		}
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->labels = 0;
		ri[i]->len = 0;
		ri[i]->mapq = -1.0;
		ri[i]->cigar = 0;
		ri[i]->errors = 0;
		ri[i]->strand = 0;
	}
	return ri;
}

void free_read_info(struct read_info** ri, int numseq)
{
	int i;
	
	for(i = 0; i < numseq;i++){
		if(ri[i]->cigar){
			MFREE(ri[i]->cigar);
		}
		
		if(ri[i]->md){
			MFREE(ri[i]->md);
		}
		
		if(ri[i]->labels){
			MFREE(ri[i]->labels);
		}
		if(ri[i]->name){
			MFREE(ri[i]->name);
		}
		
		if(ri[i]->seq){
			MFREE(ri[i]->seq);
		}
		if(ri[i]->qual){
			MFREE(ri[i]->qual );
		}
		
		MFREE(ri[i]);
	}
	MFREE(ri);
}









