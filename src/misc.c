
#include "samstat.h"
#include "nuc_code.h"

#include "misc.h"
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>

static unsigned long next = 1;


/** \var float logsum_lookup
    \brief Lookup table.
*/

/** \fn int byg_count(char* pattern,char*text)
    \brief Counts occurance of pattern text. 
    \param pattern string containing pattern.
    \param text string containing text.
    \exception Input strings must be 0 terminated. 
 
    \return number of occurences.
*/
int byg_count(char* pattern,char*text)
{
        int Tc;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }
	
        int m = (int) strlen(pattern);
        int n = (int) strlen(text);
        int count = 0;
	
        if(m > n){
                return -1;
        }
        int mb = (1 << (m-1));
	
        for (i= 0;i < m;i++){
                T[(int)toupper(pattern[i])] |= (1 << i);
        }
	
        for (i = 0;i < n;i++){
                s <<= 1;
                s |= 1;
                Tc = T[(int)toupper(text[i])];
                s &= Tc;
                if(s & mb){
                        count++;
                }
        }
        return count;
}


/** \fn int byg_end(const char* pattern,const char*text)
    \brief Finds pattern in text and returns the end coordinate of match.
    \param pattern string containing pattern.
    \param text string containing text.
    \exception Input strings must be 0 terminated.
    \return index.
*/

int byg_end(const char* pattern,const char*text)
{
        const char* tmp = 0;
        int Tc;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }
	
        int m = (int)strlen(pattern);
        int n = (int)strlen(text);
        if (m > n){
                i = m;
                m = n;
                n = i;
                tmp = text;
                text = pattern;
                pattern = tmp;
        }
	
        int mb = (1 << (m-1));
	
        for (i= 0;i < m;i++){
                T[(int)pattern[i]] |= (1 << i);
        }
	
        for (i = 0;i < n;i++){
                s <<= 1;
                s |= 1;
                if(!text[i]){
                        return -1;
                }
                Tc = T[(int)text[i]];
                s &= Tc;
                if(s & mb){
                        return i+1;
                }
        }
        return 0;
}

/** \fn int count_string(const char*p,const char** suffix,int h,int len)
    \brief Counts occurance of p in a suffix array.
    \param p string containing pattern.
    \param suffix suffix array. 
    \param h size of suffix array.
    \param len length of pattern.
    \return index.
*/

int count_string(const char*p,const char** suffix,int h,int len)
{
        int a,b;
        //for(i = 0; i < 1000000;i++){
        a = binsearch_down(p,suffix,h,len);
        b = binsearch_up(p,suffix,h,len);
        return b-a;
}

/** \fn int binsearch_down(const char*p,const char** suffix,int h,int len)
    \brief finds first occurance of p in suffix array.
    \param p string containing pattern.
    \param suffix suffix array.
    \param h size of suffix array.
    \param len length of pattern.
    \return index.
*/
int binsearch_down(const char*p,const char** suffix,int h,int len)
{
        int m = 0;
        int l = 0;
        /*if (t_long_strncmp(p,text+suffix[l],len)<= 0){
          l = l;
          }else */
        if(strncmp(p,suffix[h],len) >  0){
                return h;
        }else{
                while(h-l > 1){
                        //m = (l+h)/2;
                        m = (l + h) >> 1;
                        if(strncmp(p,suffix[m],len) <= 0){
                                h = m;
                        }else{
                                l = m;
                        }
                }
        }
        return l+1;
}

/** \fn int binsearch_up(const char*p,const char** suffix,int h,int len)
    \brief finds last occurance of p in suffix array.
    \param p string containing pattern.
    \param suffix suffix array.
    \param h size of suffix array.
    \param len length of pattern.
    \return index.
*/
int binsearch_up(const char*p,const char** suffix,int h,int len)
{
        int m = 0;
        int l = 0;
        /*if (t_long_strncmp(p,text+suffix[l],len)<= 0){
          l = l;
          }else*/
        if(strncmp(p,suffix[h],len) >  0){
                return h;
        }else{
                while(h-l > 1){
                        //m = (l+h)/2;
                        m = (l + h) >> 1;
                        if(strncmp(p,suffix[m],len) < 0){
                                h = m;
                        }else{
                                l = m;
                        }
                }
        }
        return l+1;
}



char* append_message(char* old_message, char* new_message)
{
        static size_t message_len = 0;
        struct tm *ptr;
        char time_string[200];
        int hour;
        //char am_or_pm;
        time_t current = time(NULL);
        ptr = localtime(&current);
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
	
        strftime(time_string, 200, "[%F %H:%M:%S]\t", ptr);
        fprintf(stderr,"%s%s",time_string,new_message );
        //%H:%M:%S.000
        //sprintf(time_string,"[%.2d-%.2d-%d %2d:%.2d%cm\t",ptr->tm_mon + 1,ptr->tm_mday, ptr->tm_year + 1900,hour,ptr->tm_min, am_or_pm );
        size_t time_len = strlen(time_string);
	
        size_t added_len = strlen(new_message);
	
        if(message_len == 0){
                MMALLOC(old_message,sizeof(char) *( time_len+added_len+1));
        }else{
                MREALLOC(old_message,sizeof(char) *( message_len + time_len+added_len + 1));

        }
	
	
        //char *concat = (char*) malloc(len1 + len2 + 1);
        memcpy(old_message+message_len, time_string, time_len+1);
	
        memcpy(old_message+message_len+time_len, new_message, added_len+1);
        message_len =strlen(old_message);

	
        return old_message;
ERROR:
        return NULL;
}


/** \fn int qsort_string_cmp(const void *a, const void *b)
    \brief Compares two strings. 
    Used to sort arrays of string using qsort.   
    \param a void pointer to first string. 
    \param b void pointer to second string. 
*/

int qsort_string_cmp(const void *a, const void *b)
{
        const char **one = (const char **)a;
        const char **two = (const char **)b;
        return strcmp(*one, *two);
}

/** \fn int qsort_flt_cmp(const void *a, const void *b)
    \brief Compares two floats.
    Used to sort arrays of floats.
    \param a void pointer to first float.
    \param b void pointer to second float.
*/
int qsort_flt_cmp(const void * a, const void * b)
{
        //const float a  = (float) *elem1;
        if(*(const float*)a > *(const float*)b)
                return -1;
        return *(const float*)a < *(const float*)b;
	
        //return (*(float*) a) - (*(float*) b);
}


/** \fn double gaussian_pdf(double x, double m,double s)
    Calculates the gaussian probability density function $P(X=x)$ for a normal distribution.
    \param x value.
    \param m mean.
    \param s standard deviation.
*/
double gaussian_pdf(double x, double m,double s)
{
        double a = (x-m) / s;
        return INV_SQRT_2PI / s *exp(-0.5 * a * a);
}


/** \fn unsigned int pop(int x)
    \brief Counts bits in x.
    \param x value.
    \return number of bits in \a x
*/
unsigned int pop(int x)
{
        unsigned int n = 0;
        while(x != 0){
                n = n +1;
                x = x &(x-1);
        }
        return n;
}

/** \fn int bpm(const  char* t,const  char* p,int n,int m)
    \brief Calculates edit distance between two strings. 
    \param t string1.
    \param p string 2.
    \param n length of t.
    \param m length of p.
    \return exit distance.
*/
int bpm(const  char* t,const  char* p,int n,int m)
{
        register unsigned long int i;//,c;
        unsigned long int diff;
        unsigned long int B[255];
        if(m > 31){
                m = 31;
        }
	
        unsigned long int k = m;
        //static int counter = 0;
        register unsigned long int VP,VN,D0,HN,HP,X;
	
        long int MASK = 0;
        //c = 0;
	
        diff = m;
	
        for(i = 0; i < 255;i++){
                B[i] = 0;
        }
	
        for(i = 0; i < m;i++){
                B[(int)(p[i] )] |= (1ul << i);
        }
	
        //c = 0;
        VP = 0xFFFFFFFFFFFFFFFFul;
        VN = 0ul;
        m--;
        MASK = 1ul << m;
	
        for(i = 0; i < n;i++){
                X = (B[(int)(t[i])  ] | VN);
                D0 = ((VP+(X&VP)) ^ VP) | X ;
                HN = VP & D0;
                HP = VN | ~(VP | D0);
                X = HP << 1ul;
                VN = X & D0;
                VP = (HN << 1ul) | ~(X | D0);
                diff += (HP & MASK) >> m;
                diff -= (HN & MASK) >> m;
                if(diff < k){
                        k = diff;
                        //if(k <= limit){
                        //	return (int)k;
                        //}
			
                }
		
        }
        return (int)k;
}



/** \fn int bpm(const  char* t,const  char* p,int n,int m)
    \brief Calculates edit distance between two strings.
    \param t string1.
    \param p string 2.
    \param n length of t.
    \param m length of p.
    \return exit distance.
*/
int bpm_global(const  char* t,const  char* p,int n,int m)
{
        register unsigned long int i;//,c;
        unsigned long int diff;
        unsigned long int B[255];
	
        int c;
        char* p1= 0;
        char* p2 = 0;
	
        MMALLOC(p1, sizeof(char) * (n+11));
        MMALLOC(p2, sizeof(char) * (m+11));
	
	
        for(i = 0; i < 5;i++){
                p1[i] = 'F';
                p2[i] = 'F';
        }
        c = 5;
        for(i = 0; i < n;i++){
                p1[c] = t[i];
                c++;
        }
        for(i = 0; i < 5;i++){
                p1[c] = 'Q';
                c++;
        }
        p1[c] = 0;
        n = c;
	
	
        c = 5;
        for(i = 0; i < m;i++){
                p2[c] = p[i];
                c++;
        }
        for(i = 0; i < 5;i++){
                p2[c] = 'Q';
                c++;
        }
        p2[c] = 0;
        m = c;
	
	
	
	
        if(m > 31){
                m = 31;
        }
	
        unsigned long int k = m;
        //static int counter = 0;
        register unsigned long int VP,VN,D0,HN,HP,X;
	
        long int MASK = 0;
        //c = 0;
	
        diff = m;
	
        for(i = 0; i < 255;i++){
                B[i] = 0;
        }
	
        for(i = 0; i < m;i++){
                B[(int)(p2[i] )] |= (1ul << i);
        }
	
        //c = 0;
        VP = 0xFFFFFFFFFFFFFFFFul;
        VN = 0ul;
        m--;
        MASK = 1ul << m;
	
        for(i = 0; i < n;i++){
                X = (B[(int)(p1[i])  ] | VN);
                D0 = ((VP+(X&VP)) ^ VP) | X ;
                HN = VP & D0;
                HP = VN | ~(VP | D0);
                X = HP << 1ul;
                VN = X & D0;
                VP = (HN << 1ul) | ~(X | D0);
                diff += (HP & MASK) >> m;
                diff -= (HN & MASK) >> m;
                if(diff < k){
                        k = diff;
                        //if(k <= limit){
                        //	return (int)k;
                        //}
			
                }
		
        }
        MFREE(p1);
        MFREE(p2);
	
        return (int)k;
ERROR:
        return -1;
}



/** \fn int bpm_check_error(const unsigned char* t,const unsigned char* p,int n,int m)
    \brief Calculates edit distance between two strings.
    \param t string1.
    \param p string 2.
    \param n length of t.
    \param m length of p.
    \return exit distance.
*/
int bpm_check_error(const unsigned char* t,const unsigned char* p,int n,int m)
{
        register unsigned long int i;//,c;
        unsigned long int diff;
        unsigned long int B[5];
        if(m > 31){
                m = 31;
        }
	
        unsigned long int k = m;
        //static int counter = 0;
        register unsigned long int VP,VN,D0,HN,HP,X;
	
        long int MASK = 0;
        //c = 0;
	
        diff = m;
	
        for(i = 0; i < 5;i++){
                B[i] = 0;
        }
	
        for(i = 0; i < m;i++){
                B[(int)(p[i] & 0x3)] |= (1ul << i);
        }
	
        //c = 0;
        VP = 0xFFFFFFFFFFFFFFFFul;
        VN = 0ul;
        m--;
        MASK = 1ul << m;
        for(i = 0; i < n;i++){
                X = (B[(int)(t[i] &0x3)  ] | VN);
                D0 = ((VP+(X&VP)) ^ VP) | X ;
                HN = VP & D0;
                HP = VN | ~(VP | D0);
                X = HP << 1ul;
                VN = X & D0;
                VP = (HN << 1ul) | ~(X | D0);
                diff += (HP & MASK) >> m;
                diff -= (HN & MASK) >> m;
                if(diff < k){
                        k = diff;
                        //fprintf(stderr,"%ld	%ld\n",i,k);
                        //if(k <= limit){
                        //	return (int)k;
                        //}
			
                }
        }
        return (int)k;
}




/** \fn int bpm_check_error_global(const unsigned char* t,const unsigned char* p,int n,int m)
    \brief Calculates edit distance between two strings.
    \param t string1.
    \param p string 2.
    \param n length of t.
    \param m length of p.
    \return exit distance.
*/
int bpm_check_error_global(const unsigned char* t,const unsigned char* p,int n,int m)
{
        register unsigned long int i;//,c;
        unsigned long int diff;
        unsigned long int B[5];
        if(m > 63){
                m = 63;
        }
	
        //unsigned long int k = m;
        //static int counter = 0;
        register unsigned long int VP,VN,D0,HN,HP,X;
	
        long int MASK = 0;
        //c = 0;
	
        diff = m;
	
        for(i = 0; i < 5;i++){
                B[i] = 0;
        }
	
        for(i = 0; i < m;i++){
                B[(int)(p[i] & 0x3)] |= (1ul << i);
        }
	
        //c = 0;
        VP = 0xFFFFFFFFFFFFFFFFul;
        VN = 0ul;
        m--;
        MASK = 1ul << m;
        for(i = 0; i < n;i++){
                X = (B[(int)(t[i] &0x3)  ] | VN);
                D0 = ((VP+(X&VP)) ^ VP) | X ;
                HN = VP & D0;
                HP = VN | ~(VP | D0);
                X = HP << 1ul;
                VN = X & D0;
                VP = (HN << 1ul) | ~(X | D0);
                diff += (HP & MASK) >> m;
                diff -= (HN & MASK) >> m;
                //if(diff < k){
                //	k = diff;
                //fprintf(stderr,"%ld	%ld\n",i,k);
                //if(k <= limit){
                //	return (int)k;
                //}
			
                //}
        }
        return (int)diff;
}


/** \fn int validate_bpm_sse(unsigned char**  query, int* query_lengths,unsigned char* t,int n,int num)
    \brief Calculates edit distance between four queries and one target. 
    \param query array of four strings.
    \param query_lengths lengths of query strings.
    \param t target sequence.
    \param n length of t.
    \param number of queries. 
*/

int validate_bpm_sse(unsigned char**  query, int* query_lengths,unsigned char* t,int n,int num)
{
        int i,j;
        int len = 0;
	
        unsigned int _MM_ALIGN16 nuc[16];
	
        unsigned int _MM_ALIGN16 lengths[4];
	
        __m128i VP,VN,D0,HN,HP,X,MASK,K,NOTONE,diff,zero,one;
        __m128i* nuc_p;
        __m128i xmm1,xmm2;
	
        for(i = 0; i < 16;i++){
                nuc[i] = 0ul;
        }
	
        for(i = 0; i < num;i++){
                len = query_lengths[i];
                if(len > 31){
                        len = 31;
                }
		
                lengths[i] = len;
		
		
		
                for(j = 0; j < len;j++){
                        nuc[((int)(query[i][j] & 0x3u) << 2) + i] |=  (1ul << j);// (unsigned long int)(len-1-j));
                }
        }
        nuc_p = (__m128i*) nuc;
        zero = _mm_set1_epi32(0);
        one = _mm_set1_epi32(1);
        diff =  _mm_load_si128 ( (__m128i*) lengths );  // _mm_set1_epi32(m);
        VP =  _mm_set1_epi32(0xFFFFFFFFu);
        VN =  _mm_set1_epi32(0);
        NOTONE =  _mm_set1_epi32(0xFFFFFFFF);
        K =  _mm_set1_epi32(0x7FFFFFFF);
        //POS = _mm_set1_epi32(0);
	
        for(i = 0; i< 4;i++){
                lengths[i]--;
                lengths[i] = 1 << lengths[i];
        }
	
        MASK =  _mm_load_si128 ( (__m128i*) lengths ); //  _mm_set1_epi32(1ul << m);
	
        for(i = 0; i < n ;i++){
                //fprintf(stderr,"%c",*t + 65);
                X = _mm_or_si128 (*(nuc_p +( (int)(*t)  & 0x3u) ) , VN);
                //X = (B[(int) *t] | VN);
                xmm1 = _mm_and_si128(X, VP);
                xmm2 = _mm_add_epi32(VP ,xmm1);
                xmm1 = _mm_xor_si128 (xmm2, VP);
                D0 = _mm_or_si128(xmm1, X);
                //D0 = ((VP+(X&VP)) ^ VP) | X ;
                HN = _mm_and_si128(VP, D0);
                //HN = VP & D0;
                xmm1 = _mm_or_si128(VP, D0);
                xmm2 = _mm_andnot_si128 (xmm1,NOTONE);
                HP = _mm_or_si128(VN, xmm2);
                //HP = VN | ~(VP | D0);
                X = _mm_slli_epi32(HP,1);
                //X = HP << 1ul;
                VN = _mm_and_si128(X, D0);
                //VN = X & D0;
                xmm1 = _mm_slli_epi32(HN,1);
                xmm2 = _mm_or_si128(X, D0);
                xmm2 = _mm_andnot_si128 (xmm2,NOTONE);
                VP = _mm_or_si128(xmm1, xmm2);
                //VP = (HN << 1ul) | ~(X | D0);
                xmm1 = _mm_and_si128(HP, MASK);
                xmm2 = _mm_cmpgt_epi32(xmm1, zero);
                diff = _mm_add_epi32(diff , _mm_and_si128( xmm2, one));
                //diff += (HP & MASK) >> m;
                xmm1 = _mm_and_si128(HN, MASK);
                xmm2 = _mm_cmpgt_epi32(xmm1, zero);
                diff = _mm_sub_epi32(diff,  _mm_and_si128( xmm2, one));
                //diff -= (HN & MASK) >> m;
                xmm1 = _mm_cmplt_epi32(diff, K);
                xmm2 = _mm_and_si128(xmm1, diff);
                K = _mm_or_si128(xmm2, _mm_andnot_si128  (xmm1,K));
                t++;
        }
	
        _mm_store_si128 ((__m128i*) query_lengths, K);
	
        return 1;
}

char* reverse_complement(char* p,int len)
{
        char* tmp = 0;
        MMALLOC(tmp,sizeof(char)*(len +2));
        int i,c;
        c = 0;
        for(i =len-1; i >= 0;i--){
                tmp[c] = rev_nuc_code[(int)p[i]];
                c++;
        }
        tmp[c] = 0;
        for(i= 0; i < len;i++){
                p[i] = tmp[i];
        }
        MFREE(tmp);
        return p;
ERROR:
        return NULL;
}

/** \fn void reverse_sequence(char* p,int len)
    \brief Reverses sequences.
 
    \param p nucleotide sequence.
    \param len length.
*/


void reverse_sequence(char* p,int len)
{
        int c, i, j;
	
        for (i = 0, j = len - 1; i < j; i++, j--)
        {
                c = p[i];
                p[i] = p[j];
                p[j] = c;
        }
}



/* RAND_MAX assumed to be 32767 */
int myrand(void)
{
        next = next * 1103515245 + 12345;
        return((unsigned)(next/65536) % 32768);
}

void mysrand(unsigned seed)
{
        next = seed;
}

int file_exists(char* name)
{
        struct stat buf;
        int ret,local_ret;
        ret = 0;
	
        local_ret= stat ( name, &buf );
        /* File found */
        if ( local_ret == 0 )
        {
                ret++;
                //return 1;
        }
        return ret;
}




