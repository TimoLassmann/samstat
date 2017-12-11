#ifndef tagdust2_misc_h
#define tagdust2_misc_h

#include <stdio.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"


#define MAX_LINE 10000

#define SQRT2M_PI 2.506628274631


#define INV_SQRT_2PI 0.3989422804014327


#ifndef _MM_ALIGN16
#ifdef __GNUC__
#define _MM_ALIGN16 __attribute__((aligned (16)))
#endif
#ifdef __MSVC__
#define _MM_ALIGN16 __declspec(align(16))
#endif
#endif


#endif

#ifdef RTEST
#define rand() myrand()
#define srand(x) mysrand(x)

#endif



int myrand(void);
void mysrand(unsigned seed);

int byg_end(const char* pattern,const char*text);

int byg_count(char* pattern,char*text);

char* append_message(char* old_message, char* new_message);
void init_logsum();
float logsum(float a,float b);

float logsum_print(float a,float b);
unsigned int pop(int x);
float prob2scaledprob(float p);
float scaledprob2prob(float p);

int qsort_string_cmp(const void *a, const void *b);
int qsort_flt_cmp(const void * a, const void * b);

int bpm(const char* t,const char* p,int n,int m);
int bpm_global(const  char* t,const  char* p,int n,int m);


int binsearch_down(const char*p,const char** suffix,int h,int len);
int binsearch_up(const char*p,const char** suffix,int h,int len);

int count_string(const char*p,const char** suffix,int h,int len);

double gaussian_pdf(double x, double m,double s);

int validate_bpm_sse(unsigned char**  query, int* query_lengths,unsigned char* t,int n,int num);
int bpm_check_error(const unsigned char* t,const unsigned char* p,int n,int m);
int bpm_check_error_global(const unsigned char* t,const unsigned char* p,int n,int m);
char* reverse_complement(char* p,int len);
void reverse_sequence(char* p,int len);
int file_exists(char* name);

