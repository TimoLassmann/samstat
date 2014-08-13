//
//  pst.c
//  tagdust2
//
//  Created by lassmann on 2/7/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include "malloc_macro.h"
#include "nuc_code.h"
#include "misc.h"
#include "io.h"
#include <time.h>
#include "pst.h"


#include <float.h>



struct pst_node* build_pst(struct pst* pst,struct pst_node* n )
{
	char alphabet[] = "ACGTN";
	
	char tmp[MAX_PST_LEN+5];
	int i;
	int j;
	int c;
	int add;
	int len = (int) strlen(n->label);
	double sum = 0.0f;
	
	
	
	double tmp_counts_s[5];
	
	
	
	//fprintf(stderr,"NODE: %s\n", n->label);
	//for(i = 0;i < 5;i++){
	//	fprintf(stderr,"%c+%s\t%f\n",alphabet[i], n->label, n->nuc_probability[i]);
	//}
	
	//step 2 test expansion
	
	//loop though letters at present node
	if(len + 1 < MAX_PST_LEN ){
		/// search for all strings and record probabilities S+ACGT...
		/// don't search rare strings...
		/// - super mega simple ...
		
		
		for(i = 0; i < 5;i++){
			if(n->nuc_probability[i] >= pst->p_min){ /// string (suffix  + letter is not rare compared to all other letters  ) // e.g. acgC = 0.2 and not acgC = 0.0000001. hence it makes sense to "extend the suffix - test [A,C,G,T] - acgC
								
				
				
				
				
				//init longer suffix
				tmp[0] = alphabet[i];
				for(j = 1; j < len+1;j++){
					tmp[j] = n->label[j-1];
				}
				
				sum = 0.0;
				for(j = 0; j < 5;j++){
					tmp[len+1]  = alphabet[j];
					tmp[len+2] = 0;
					c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,len+2);
					tmp_counts_s[j] = c;
					sum+= c;
				}
				
				add = 0;
				for(j = 0; j < 5;j++){
					if(tmp_counts_s[j] /  (pst->numseq * (pst->mean_length - (len+1))) > 0.001){
						add = 1;
						break;
					}
				}
				if(add){
					// here I know that probablility of 'X' is non-neglible AND that there existsa string 'X' - [A,C,G,T,N] which is frequent in the data - hence I add...
					n->next[i] = alloc_node(n->next[i] ,tmp,len+1);
					add = 0;
					for(j = 0; j < 5;j++){
						if((tmp_counts_s[j]/sum) / n->nuc_probability[j] >= pst->r){
							add++;
						}
						
						if((tmp_counts_s[j]/sum) / n->nuc_probability[j] <= 1.0/ pst->r){
							add++;
						}
						
						n->next[i]->nuc_probability[j] = tmp_counts_s[j]/sum;
							
					}
					if(add){
						n->next[i]->in_T = 1;
					}
					n->next[i] = build_pst(pst,n->next[i]  );

				}
			}
		}
	}
	c= 0;
	for(i = 0; i < 5;i++){
		if(n->next[i]){
			c+= n->next[i]->in_T;
		}
	}
	if(c){
		n->in_T = 1;
	}
	return n;
}



struct pst_node* build_ppt(struct pst* pst,struct pst_node* n )
{
	char alphabet[] = "ACGTN";
	
	char tmp[MAX_PST_LEN+5];
	int i;
	int j;
	int c;
	int add;
	int len = (int) strlen(n->label);
	double sum = 0.0f;
	
	
	
	double tmp_counts_s[5];
	
	
	
	//fprintf(stderr,"NODE: %s\n", n->label);
	//for(i = 0;i < 5;i++){
	//	fprintf(stderr,"%c+%s\t%f\n",alphabet[i], n->label, n->nuc_probability[i]);
	//}
	
	//step 2 test expansion
	
	//loop though letters at present node
	if(len + 1 < MAX_PST_LEN ){
		/// search for all strings and record probabilities S+ACGT...
		/// don't search rare strings...
		/// - super mega simple ...
		
		
		for(i = 0; i < 5;i++){
			if(n->nuc_probability[i] >= pst->p_min){ /// string (suffix  + letter is not rare compared to all other letters  ) // e.g. acgC = 0.2 and not acgC = 0.0000001. hence it makes sense to "extend the suffix - test [A,C,G,T] - acgC
				
				
				
				
				
				//init longer prefix!!!! 
				
				for(j = 0; j < len;j++){
					tmp[j+1] = n->label[j];
				}
				tmp[len+1] = alphabet[i];
				
				sum = 0.0;
				for(j = 0; j < 5;j++){
					tmp[0]  = alphabet[j];
					tmp[len+2] = 0;
					c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,len+2);
					tmp_counts_s[j] = c;
					sum+= c;
				}
				
				add = 0;
				for(j = 0; j < 5;j++){
					if(tmp_counts_s[j] /  (pst->numseq * (pst->mean_length - (len+1))) > 0.001){
						add = 1;
						break;
					}
				}
				if(add){
					// here I know that probablility of 'X' is non-neglible AND that there existsa string 'X' - [A,C,G,T,N] which is frequent in the data - hence I add...
					n->next[i] = alloc_node(n->next[i],tmp+1,len+1);
					add = 0;
					for(j = 0; j < 5;j++){
						if((tmp_counts_s[j]/sum) / n->nuc_probability[j] >= pst->r){
							add++;
						}
						
						if((tmp_counts_s[j]/sum) / n->nuc_probability[j] <= 1.0/ pst->r){
							add++;
						}
						
						n->next[i]->nuc_probability[j] = tmp_counts_s[j]/sum;
						
					}
					if(add){
						n->next[i]->in_T = 1;
					}
					n->next[i] = build_ppt(pst,n->next[i]);
					
				}
			}
		}
	}
	c= 0;
	for(i = 0; i < 5;i++){
		if(n->next[i]){
			c+= n->next[i]->in_T;
		}
	}
	if(c){
		n->in_T = 1;
	}
	
	return n;
}

struct pst* alloc_pst(int numseq)
{
	struct pst* pst = NULL;
	
	MMALLOC(pst, sizeof(struct pst));
	
	pst->suffix_array = NULL;
	pst->L = MAX_PST_LEN;
	pst->alpha = 0.0f;
	pst->p_min = 0.0001f;
	pst->lamba = 0.001f;
	pst->r = 1.05f;
	pst->total_len = 0;
	pst->pst_root = NULL;
	pst->ppt_root = NULL;
	pst->rank_array = 0;

	
	
	
	pst->current_suffix_size =numseq;
	MMALLOC(pst->suffix_array,sizeof(char*)* pst->current_suffix_size);
	
	pst->pst_root = alloc_node(pst->pst_root,"",0);
	pst->ppt_root = alloc_node(pst->ppt_root,"",0);
	

	return pst;
}

struct pst_node* alloc_node(struct pst_node* n,char* string,int len)
{
	int i;
	MMALLOC(n ,sizeof(struct pst_node));
	
	//assert(n!=0);
	
	n->label = malloc(sizeof(char) *(len+1));
	//assert(n->label != 0);
	
	for(i = 0; i < len;i++){
		n->label[i] =string[i];
	}
	n->label[len] = 0;
	n->in_T = 0;
	n->occ = 0;
	n->last_seen = -1;
	//n->bit_occ = 0;
	for(i =0; i < 5;i++){
		n->next[i] = 0;
		n->nuc_probability[i] = 0.2f;
	}

	return n;
}


struct read_info**  scan_read_with_pst(struct read_info** ri,struct pst* pst)
{
	int i,j;
	float P_T;
	//float P_PT;
	float P_R = 0;
	
	float* base_p = pst->pst_root->nuc_probability;
	char* seq;
	
	float total_T = prob2scaledprob(1.0);
	float total_R = prob2scaledprob(1.0);
	
	float A,B;
	
	//scan to cpount occurances... 
	for(i = 0; i < pst->numseq;i++){
		
		seq = ri[i]->seq;
		for(j = 0; j < ri[i]->len; j++ ){
			pst->pst_root = count_pst_lables(pst->pst_root, seq,  j, i);
			pst->ppt_root = count_ppt_lables(pst->ppt_root, seq, j , i);
			
		}
	}
	
	//fprintf(stderr,"%f\n",pst->numseq);
	
	for(i = 0; i < pst->numseq;i++){
		if(!ri[i]->qual){
			ri[i]->qual = malloc(sizeof(char)* (ri[i]->len+1));
		}
		for(j = 0; j < ri[i]->len; j++ ){
			ri[i]->qual[j] = 48;
		}
		//qual = ri[i]->qual;
		seq = ri[i]->seq;
		
		//if(!i){
		
		//}
		
		P_T = prob2scaledprob(1.0);
		P_R = prob2scaledprob(1.0);
		//P_PT =  prob2scaledprob(1.0);
		//fprintf(stdout,"%s\n",seq );
				
		for(j = 0; j < ri[i]->len; j++ ){
			P_R = P_R + prob2scaledprob(base_p[nuc_code[(int)seq[j]]]);
			
			
			//get_occ
			A = get_pst_prob(pst->pst_root, seq,  nuc_code[(int)seq[j]], j, i);
			B = get_ppt_prob(pst->ppt_root, seq,  nuc_code[(int)seq[j]], j, i);
			
			
			P_T = P_T + prob2scaledprob( A > B ? A:B );
			
			//fprintf(stdout,"%d	%f	%f\n", j , A,B);
			ri[i]->qual[j] = 48 + (int) (10.0* ( A > B ? A:B ));
		}
		
		total_T = total_T + P_T;
		total_R = total_R + P_R;
		
		ri[i]->mapq = P_T-P_R;
	}
	return ri;
}



int add_patterns(struct pst_node** all_patterns, struct pst_node* n,int num)
{
	int i;
	
	int internal = 0;
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			internal++;
		}
	}
	if (!internal){
	if(strlen(n->label) > 6){
		if(n->in_T){
			all_patterns[num] = n;
			num += 1;
			//fprintf(stderr,"ADDED: %s	%d\n",n->label,num);
		}
	}
		
	}
	
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			//fprintf(stderr,"Going:%d\n",i);
			num = add_patterns(all_patterns,n->next[i],num);
		}
	}
	
	return num;
}

int count_patterns(struct pst_node* n,int num)
{
	int i;
	if(strlen(n->label) > 2){
		if(n->in_T){
			num += 1;
			//fprintf(stderr,"ADDED: %s	%d\n",n->label,num);
		}
	}
	
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			//fprintf(stderr,"Going:%d\n",i);
			num = count_patterns(n->next[i],num);
		}
	}

	return num;
}

float get_pst_prob(struct pst_node* n, char* string,int target, int pos,int seq_id)
{
	int c;
	//if(!seq_id){
	//fprintf(stderr,"%d	%s	%f	(%d)	- OCC:%d\n", target,n->label,n->nuc_probability[target],n->in_T, n->occ);
	//}
	/*if(n->last_seen != seq_id){
		n->occ++;
		n->last_seen = seq_id;
	}*/
	if(pos == 0){
		
		return n->nuc_probability[target];
	}
	pos = pos -1;
	c = nuc_code[(int)string[pos]];
	if(n->next[c]){
		return get_pst_prob(n->next[c], string, target,pos,seq_id);
	}else{
		
		return n->nuc_probability[target];
	}
}




float get_ppt_prob(struct pst_node* n, char* string,int target, int pos,int seq_id)
{
	int c;
	//if(!seq_id){
	//fprintf(stderr,"%d	%s	%f	(%d)	- OCC:%d\n", target,n->label,n->nuc_probability[target],n->in_T, n->occ);
	//}
	/*if(n->last_seen != seq_id){
	 n->occ++;
	 n->last_seen = seq_id;
	 }*/
	if(string[pos+1] == 0){
		
		return n->nuc_probability[target];
	}
	pos = pos +1;
	c = nuc_code[(int)string[pos]];
	if(n->next[c]){
		return get_ppt_prob(n->next[c], string, target,pos,seq_id);
	}else{
		
		return n->nuc_probability[target];
	}
}

int get_occ(struct pst_node* n, char* string,int target, int pos,int seq_id)
{
	//if(!seq_id){
	//fprintf(stderr,"%d	%s	%f	(%d)	- OCC:%d\n", target,n->label,n->nuc_probability[target],n->in_T, n->occ);
	//}
	/*if(n->last_seen != seq_id){
	 n->occ++;
	 n->last_seen = seq_id;
	 }*/
	if(pos == 0){
		return n->occ;
	}
	pos = pos -1;
	int c = nuc_code[(int)string[pos]];
	if(n->next[c]){
		return get_occ(n->next[c], string, target,pos,seq_id);
	}else{
		return n->occ;
	}
}


struct pst_node*  count_pst_lables(struct pst_node* n, char* string, int pos,int seq_id)
{
	//if(!seq_id){
	//fprintf(stderr,"%d	%s	%f	(%d)	- OCC:%d\n", target,n->label,n->nuc_probability[target],n->in_T, n->occ);
	//}
	
	int c = nuc_code[(int)string[pos]];
	if(n->next[c]){
		if(n->next[c]->last_seen != seq_id){
			if(!bit_test(n->next[c]->bit_occ , seq_id)){
				n->next[c]->occ++;
				n->next[c]->bit_occ = bit_set(n->next[c]->bit_occ , seq_id);
			}
			n->next[c]->last_seen = seq_id;
		}
		pos = pos -1;
		if(pos != -1){
			n->next[c] =  count_pst_lables(n->next[c], string,pos,seq_id);
		}
	}else{
		return n;
	}
	
	/*if(n->last_seen != seq_id){
		n->occ++;
		//n->bit_occ = bit_set(n->bit_occ, seq_id);
		n->last_seen = seq_id;
	}
	//fprintf(stderr,"%s	(L-PST)\n", n->label );
	//fprintf(stderr,"%s\n", string+pos);
	//fprintf(stderr,"%s\n", string);
	
	if(pos == 0){
		return n;
	}
	pos = pos -1;
	
	if(n->next[c]){
		n->next[c] =  count_pst_lables(n->next[c], string,pos,seq_id);
	}else{
		return n;
	}*/
	return n;
}


struct pst_node*  count_ppt_lables(struct pst_node* n, char* string, int pos,int seq_id)
{
	int c = nuc_code[(int)string[pos]];
	if(!string[pos]){
		return n;
	}
	if(n->next[c]){
		if(n->next[c]->last_seen != seq_id){
			if(!bit_test(n->next[c]->bit_occ , seq_id)){
				n->next[c]->occ++;
				n->next[c]->bit_occ = bit_set(n->next[c]->bit_occ , seq_id);
			}
			n->next[c]->last_seen = seq_id;
		}
		pos = pos +1;
		//n = n->next;
		//if(pos){
			n->next[c] =  count_ppt_lables(n->next[c], string,pos,seq_id);
		//}
	}else{
		return n;
	}

	
	return n;
}




void free_pst(struct pst_node* n)
{
	int i;
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			free_pst(n->next[i]);
		}
	}
	if(n->bit_occ){
		free(n->bit_occ);
	}
	free(n->label);
	free(n);

}

void print_pst(struct pst* pst,struct pst_node* n, struct read_info** ri )
{
	int i;
	int internal;
	
	double p ,e;
	
//	int len = (int)strlen(n->label);
	
//	double N1,N2,U1,U2,R1,R2,Z;
//	int c;
//	struct ranks** rank_array  = pst->rank_array;
	//struct ranks** rank_array = malloc(sizeof(struct ranks*) * (int)(pst->numseq));
	//for(i = 0; i  <pst->numseq;i++){
	//	rank_array[i] = malloc(sizeof(struct ranks));
	//}
	//char alphabet[] = "ACGTN";
	//if(strlen(n->label) > 2){
		internal = 0;
		for(i = 0;i < 5;i++){
			if(n->next[i]){
				internal++;
			}
		}
		if(!internal){
			/*
			//fprintf(stderr,"%p\n",n);
			//fprintf(stderr,"%s	%d	%d\n", n->label,n->in_T, count_string(n->label,(const char**)pst->suffix_array,pst->suffix_len-1,len));
			c = 0;
			N1 = 0;
			N2 = 0;
			//fprintf(stderr,"NUMSEQ:::::%f\n",pst->numseq );
			for(i = 0 ;i  < pst->numseq;i++){
				//fprintf(stderr,"%d ",i);
				if(bit_test(n->bit_occ, i)){
					rank_array[c]->sample = 0;
					rank_array[c]->value = ri[i]->mapq;
					N1++;
				}else{
					rank_array[c]->sample = 1;
					rank_array[c]->value = ri[i]->mapq;
					N2++;
				}
				c++;
				
			}
			
			
			qsort((void *)  rank_array, N1+N2, sizeof(struct ranks* ),(compfn) establish_rank);
			
			U1 = 0.0;
			R1 = 0.0;
			
			U2 = 0.0;
			R2 = 0.0;
			for(i = 0;i < N1+N2;i++){
				
				if(!rank_array[i]->sample){
					R1 += (i+1);
				}else{
					R2 += (i+1);
				}
				//if(i < 10){
				//	fprintf(stderr,"%d\t%d\t%f\n",i, rank_array[i]->sample,rank_array[i]->value);
				//}
			}
			
			U1 = R1 - (N1*(N1 + 1.0))/2.0;
			U2 = R2 - (N2*(N2 + 1.0))/2.0;
			
			if(U2 < U1){
				U1 = U2;
			}
			
			
			Z = (U1- (N1*N2 /2.0) )/  sqrt((N1 * N2 *(N1 + N2 +1) )/  12.0 );
			//	fprintf(stderr,"%f	%f	%f	%f	\n",N1,N2,U1,U2);
			//Z = 1;
			N1 = prob2scaledprob(1.0);
			for(i = 0 ;i < strlen(n->label);i++){
				N1 += prob2scaledprob(pst->pst_root->nuc_probability[nuc_code5[ (int)n->label[i]]] );
			}
			N1 = N1 + prob2scaledprob( pst->mean_length - strlen(n->label)) + prob2scaledprob(pst->numseq) ;
			*/
			
		p =  (double) n->occ / (double)  pst->numseq;
		
		e = p * log2(p);
		
		p = 1.0 -  (double) n->occ / (double)  pst->numseq;
		e+=  p * log2(p);
		e *= -1.0;
			fprintf(stderr,"%s	%d	%d	Entropy:%f	%f\n", n->label,n->in_T,n->occ,  e,pst->numseq);
			//for(i = 0;i < 5;i++){
				//if(n->next[i]){
			//	fprintf(stderr,"%f ",n->nuc_probability[i]);
				//}else{
				//	fprintf(stderr,"%f S\t",n->nuc_probability[i]);
				//}
			//}
			
			
		//	fprintf(stderr,"\n");
		}
	//}
	
	
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			if(n->next[i]->in_T){
			//fprintf(stderr,"Going:%d\n",i);
			print_pst(pst,n->next[i],ri);
			}
		}
	}
}


struct pst_node* alloc_bit_occ_pst(struct pst_node* n, int num)
{
	int i;
	n->bit_occ = _mm_malloc(sizeof(int)* (1+ num / BITSPERWORD) , 16);
	
	//assert(n->bit_occ != 0);
	
	for(i = 0; i < (1+ num / BITSPERWORD);i++){
		n->bit_occ[i] = 0;
	}
	
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			n->next[i] = alloc_bit_occ_pst(n->next[i],num);
		}
	}
	return n;
}


int* bit_set(int*a, int i)
{
	a[i >> SHIFT] |= (1 <<(i & MASK));
	return a;
}


int* bit_clr(int*a, int i)
{
	a[i >> SHIFT] &= ~(1 <<(i & MASK));
	return a;
}

int bit_test(int*a, int i)
{
	return a[i >> SHIFT] & (1 <<(i & MASK));
}


int establish_rank(const void *a, const void *b)
{
	struct ranks* const *one = a;
	struct ranks* const *two = b;
	
	if((*one)->value >  (*two)->value){
		return -1;
	}else{
		return 1;
	}
}


int sort_pst_nodel_according_to_label(const void *a, const void *b)
{
	struct pst_node * const *one = a;
	struct pst_node* const *two = b;
	
	return strcmp((*one)->label,(*two)->label) ;
}



int sort_pst_nodel_according_to_occ(const void *a, const void *b)
{
	struct pst_node * const *one = a;
	struct pst_node* const *two = b;
	
	if((*one)->occ >  (*two)->occ){
		return -1;
	}else{
		return 1;
	}
	//return strcmp((*one)->label,(*two)->label) ;
}








int qsort_suffix_node_string_cmp(const void *a, const void *b)
{
	struct suffix_node* const *one  = a;
	struct suffix_node* const *two  = b;
	
	return strcmp((*one)->string,(*two)->string) ;
}






