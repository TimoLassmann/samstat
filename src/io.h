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

#define SEEK_START 0
#define SEEK_END 2

#include <unistd.h>

#include "interface.h"

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



