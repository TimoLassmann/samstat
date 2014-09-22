#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(SAMstat,  1.5, timolassmann@gmail.com)
AC_CONFIG_HEADER([src/config.h])
AM_INIT_AUTOMAKE([ -Wno-portability ])

CXXFLAGS="$CXXFLAGS"

AC_ARG_ENABLE(debug,
              [  --enable-debug          Enable debugging information],
              USE_DEBUG="$enableval", USE_DEBUG="no")

if test $USE_DEBUG = yes ; then 
	DEBUG=1
	CFLAGS="$CFLAGS -D DEBUG -ggdb -Wall -m64 -std=gnu99 "
else
	CFLAGS="$CFLAGS -O3 -Wall -m64 -std=gnu99"
	DEBUG=0
fi


AC_ARG_ENABLE(valgrind-tests,
    [--enable-valgrind-tests run self tests under valgrind],
    opt_valgrind_tests="$enableval", opt_valgrind_tests="no")
    

OPTS="-q --trace-children=yes --error-exitcode=1 --dsymutil=yes --leak-check=yes"
  # Run self-tests under valgrind?
if test "$opt_valgrind_tests" = "yes"; then
	CFLAGS="-DDEBUG -ggdb -pedantic -Wall -m64 -std=gnu99"
	AC_CHECK_PROGS(VALGRIND, valgrind)
	VALGRIND="$VALGRIND $OPTS"
else
	VALGRIND=
fi

AC_MSG_CHECKING([whether self tests are run under valgrind])
AC_MSG_RESULT($opt_valgrind_tests)

AC_CHECK_PROGS([SAMTOOLS], [samtools])
if test -z "$SAMTOOLS"; then
	AC_MSG_ERROR([

********************************************************
* The samtools package is not installed on your system.
* Download samtools from here: 
* http://samtools.sourceforge.net
* or 
* check you PATH variable.
********************************************************

])
fi


# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([ ctype.h getopt.h  math.h  stdio.h stdlib.h string.h sys/stat.h time.h])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])

AC_OUTPUT


echo "
     Configuration:
     
     compiler: ${CFLAGS}
     valgrind: ${VALGRIND}

"



