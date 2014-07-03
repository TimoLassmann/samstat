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

/*! \file nuc_code.c
 \brief Initializes nucleotide conversion arrays. 
 
 Initializes nucleotide arrays.
 \author Timo Lassmann
 \bug No known bugs.
 */
#include <stdlib.h>
#include "nuc_code.h"


/** \fn void init_nuc_code()
 \brief Initializes nucleotide conversion arrays. 
 
 */
void init_nuc_code()
{
	int i;
	for(i = 0;i < 256;i++){
		nuc_code[i] = 4;
	}
	
	
	nuc_code[46] = 5;/// dot - will initialize N to p = 1
	
	nuc_code[65] = 0;//A Adenine
	nuc_code[67] = 1;//C	Cytosine
	nuc_code[71] = 2;//G	Guanine
	nuc_code[84] = 3;//T	Thymine
	nuc_code[85] = 3;//U	Uracil

	nuc_code[65+32] = 0;//A Adenine
	nuc_code[67+32] = 1;//C	Cytosine
	nuc_code[71+32] = 2;//G	Guanine
	nuc_code[84+32] = 3;//T	Thymine
	nuc_code[85+32] = 3;//U	Uracil
	
	rev_nuc_code[0] = 3;//A Adenine
	rev_nuc_code[1] = 2;//C	Cytosine
	rev_nuc_code[2] = 1;//G	Guanine
	rev_nuc_code[3] = 0;//T	Thymine
	rev_nuc_code[4] = 4;//U	Uracil
}


