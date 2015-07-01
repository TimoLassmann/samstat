
// These are macros to catch errors
// All functions are adopted from the easel library (part of hmmer).

#ifndef kslib_included

#define kslib_included

#include <stdlib.h>
#include <stdio.h>		/* for FILE */
#include <stdarg.h>

#define kslOK              0    /* no error/success             */
#define kslFAIL            1    /* failure                      */
#define kslEMEM         2 /* memory failure */
#define kslEWRT         3 /* write failure */
#define kslETRD         3 /* thread failure */


#define MMALLOC(p, size) do {\
if (p != NULL){\
status =kslEMEM; \
ksl_exception(kslEMEM, FALSE, __FILE__, __LINE__, "malloc on a nun-null pointer"); \
goto ERROR;\
}\
if (((p) = malloc(size)) == NULL && (size)) {	\
status = kslEMEM;\
ksl_exception(kslEMEM, FALSE, __FILE__, __LINE__, "malloc of size %d failed", size); \
goto ERROR;\
}} while (0)


#define MREALLOC(p, newsize) do {\
void *tmpp;\
if ((p) == NULL) { (tmpp) = malloc(newsize);         }\
else             { (tmpp) = realloc((p), (newsize)); }\
if ((tmpp) != NULL) (p) = (tmpp);\
else {\
status = kslEMEM;\
ksl_exception(kslEMEM, FALSE, __FILE__, __LINE__, "realloc for size %d failed", newsize); \
goto ERROR;\
}} while (0)

#define MFREE(p) do {\
free(p);\
p = NULL;\
} while (0)


#define kslibERRBUFSIZE 128

#define kslibMSGBUFSIZE 1024

#define KSLIB_FAIL(code, errbuf, ...) do {				\
if (errbuf != NULL) snprintf(errbuf, kslibERRBUFSIZE, __VA_ARGS__);	\
return code; }							\
while (0)

#define KSLIB_XFAIL(code, errbuf, ...) do {				\
status = code;							\
if (errbuf != NULL) snprintf(errbuf, kslibERRBUFSIZE, __VA_ARGS__);	\
goto ERROR; }							\
while (0)

#define KSLIB_EXCEPTION(code, ...)  do {					\
ksl_exception(code, FALSE, __FILE__, __LINE__, __VA_ARGS__);	\
return code; }							\
while (0)

#define KSLIB_XEXCEPTION(code, ...)  do {					\
status = code;							\
ksl_exception(code, FALSE, __FILE__, __LINE__, __VA_ARGS__);	\
goto ERROR; }							\
while (0)

#define KSLIB_EXCEPTION_SYS(code, ...) do {				\
ksl_exception(code, TRUE, __FILE__, __LINE__, __VA_ARGS__);		\
return code; }							\
while (0)

#define KSLIB_XEXCEPTION_SYS(code, ...)  do {				\
status = code;							\
ksl_exception(code, TRUE, __FILE__, __LINE__, __VA_ARGS__);	\
goto ERROR; }							\
while (0)


#define KSLIB_MESSAGE(code, ...)  do {					\
ksl_exception(code, FALSE, __FILE__, __LINE__, __VA_ARGS__);	\
}							\
while (0)


#define KSLIB_snprintf(target, size, source)

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif


#if kslDEBUGLEVEL >= 1		/* for KSL_DASSERT() macros */
#include <assert.h>
#endif

#if (kslDEBUGLEVEL >= 1)
#define KSL_DPRINTF1(x)  printf x
#define KSL_DASSERT1(x)  assert x
#else
#define KSL_DPRINTF1(x)
#define KSL_DASSERT1(x)
#endif
#if (kslDEBUGLEVEL >= 2)
#define KSL_DPRINTF2(x)  printf x
#define KSL_DASSERT2(x)  assert x
#else
#define KSL_DPRINTF2(x)
#define KSL_DASSERT2(x)
#endif
#if (kslDEBUGLEVEL >= 3)
#define KSL_DPRINTF3(x)  printf x
#define KSL_DASSERT3(x)  assert x
#else
#define KSL_DPRINTF3(x)
#define KSL_DASSERT3(x)
#endif


#define MACRO_MIN(a,b)          (((a)<(b))?(a):(b))
#define MACRO_MAX(a,b)          (((a)>(b))?(a):(b))

extern void ksl_exception(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, ...);
extern void ksl_message(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, ...);
extern void ksl_error(const char *format, ...);

#endif


