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

/*! \file interface.h 
 \brief 
 
 Initializes nucleotide alphabet needed to parse input. Calls parameter parser. Calls functions to process the data. \author Timo Lassmann \bug No known bugs.
 */


#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>


struct parameters {
	char** infile; /**< @brief Names of input files. */
	char* outfile;
	int infiles;/**<  @brief Number of input files. */
	int quiet_flag;
	int num_query;/**< @brief Number of sequences to read at one time. */
	char* format;
	char* filter;
	char* train;
	char* exact5;
	char* messages;
	char* buffer;
	int gzipped;
	int bzipped;
	int dust;
	int sam;
	int fasta;
};

struct parameters* interface(struct parameters* param,int argc, char *argv[]);
void free_param(struct parameters* param);
void usage(void);



