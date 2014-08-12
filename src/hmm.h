
#ifndef MMALLOC
#include "malloc_macro.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include "math.h"

#define STARTSTATE 0
#define ENDSTATE 1


#define LOGSUM_SIZE 160000

#define SCALE 10000.0

#define MIN(a,b)          (((a)<(b))?(a):(b))
#define MAX(a,b)          (((a)>(b))?(a):(b))

#define MODE_BAUM_WELCH 0
#define MODE_FORWARD 1




struct hmm{
	float** transitions;
	float** emissions;
	float** transitions_e;
	float** emissions_e;
	float** F;
	float** B;
	int** tindex;
	
	float* F_memory;
	float* B_memory;
	
	void* data;
	
	float f_score;
	float b_score;
	int num_states;
	int alphabet_len;
	int max_seq_len;
	
};


struct hmm_data{
	char** string;
	int* length;
	int num_seq;
	float* score;
	float* weight;
	int num_threads;
	int run_mode;
	int iterations;
};


struct thread_data{
	struct hmm* hmm;
	struct parameters* param;
	struct hmm_data* data;
	int numseq;
	int start;
	int end;
};


/* convenience fiunctions */
struct hmm* run_EM_iterations (struct hmm* hmm,struct hmm_data* data);

/* Main driver functions */

struct hmm* run_pHMM(struct hmm* hmm,struct hmm_data* data);
void* do_baum_welch(void *threadarg);
void* do_forward(void *threadarg);

/* core functions */
struct hmm* backward(struct hmm* hmm, char* a, int len);
struct hmm* forward(struct hmm* hmm, char* a, int len);
struct hmm* collect_estimated(struct hmm* hmm, char* a,float weight, int len);
struct hmm* reestimate_hmm_parameters(struct hmm* hmm);


/* hmm manipulations */
struct hmm* copy_hmm(struct hmm* org);
struct hmm* copy_estimated_parameters(struct hmm* target,struct hmm* source );

/* printing */
void print_dyn_matrix(struct hmm* hmm, int seq_len);
void print_hmm_parameters(struct hmm* hmm);
void print_hmm_estimated_parameters(struct hmm* hmm);

/* mallocing */
struct hmm* malloc_hmm(int num_states, int alphabet_len, int max_seq_len);
void free_hmm(struct hmm* hmm);




struct hmm* init_hmm_simple_ID(struct hmm* hmm);




/* misc math functions */

void init_logsum();

float logsum(const float a,const float b);
//float logsub(const float a,const float b);
float prob2scaledprob(float p);
float scaledprob2prob(float p);
//float logsub_error(float a, float b);
//float logsum_error(float a, float b);

#ifdef ITEST
#include <time.h>

struct hmm* complicated_hmm(void);
void print_max_posterior(struct hmm* hmm, char* a, int len);
#endif






