
// These are macros to catch some additional errors
// All fiunctions are adopted feom the easel library (part of hmmer).

#define MMALLOC(p, size) do {\
if(p != NULL){\
fprintf(stderr, "%s:%d:%s(): pointer to be malloced is not null \n", __FILE__, \
__LINE__, __func__);\
exit(EXIT_FAILURE); \
}\
if (((p) = malloc(size)) == NULL) {\
fprintf(stderr, "%s:%d:%s(): MALLOC failed\n", __FILE__, \
__LINE__, __func__);\
exit(EXIT_FAILURE); \
}} while (0)



#define MREALLOC(p, tmp, newsize) do {\
if ((p) == NULL) {\
	(tmp) = malloc(newsize); \
} else {\
	(tmp) = realloc((p), (newsize));\
}\
if ((tmp) != NULL){\
	(p) = (tmp);\
} else {\
fprintf(stderr, "%s:%d:%s(): REALLOC failed\n", __FILE__, \
__LINE__, __func__);\
exit(EXIT_FAILURE); \
}} while (0)


#define MFREE(p) do {\
free(p);\
p = NULL;\
} while (0)


