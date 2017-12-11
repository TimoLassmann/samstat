#include <stdlib.h>
#include "nuc_code.h"

int init_nuc_code()
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
        return OK;
}


