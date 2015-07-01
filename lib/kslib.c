

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>


void ksl_exception(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, ...)
{
	va_list argp;	

	fprintf(stderr, "Fatal exception (source file %s, line %d):\n", sourcefile, sourceline);
	va_start(argp, format);
	vfprintf(stderr, format, argp);
	va_end(argp);
	fprintf(stderr, "\n");
	if (use_errno && errno) perror("system error");
	fflush(stderr);

	abort();
}


void ksl_message(int errcode, int use_errno, char *sourcefile, int sourceline, char *format, ...)
{
	va_list argp;
	
	fprintf(stderr, "Fatal exception (source file %s, line %d):\n", sourcefile, sourceline);
	va_start(argp, format);
	vfprintf(stderr, format, argp);
	va_end(argp);
	fprintf(stderr, "\n");
	if (use_errno && errno) perror("system error");
	
	fflush(stderr);

}


void ksl_error(const char *format, ...)
{
	fprintf(stderr, "\nError: ");
	va_list argp;
	va_start(argp, format);
	vfprintf(stderr, format, argp);
	va_end(argp);
	fprintf(stderr, "\n");
	fflush(stderr);
	exit(EXIT_FAILURE);
}



