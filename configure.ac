#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(SAMstat,  1.5, timolassmann@gmail.com)
AM_INIT_AUTOMAKE([ -Wno-portability ])

CXXFLAGS="$CXXFLAGS"

AC_ARG_ENABLE(debugging,
[AS_HELP_STRING([--enable-debugging],[include debugging code])
AS_HELP_STRING([--enable-debugging=x],[also set diagnostics verbosity level to <x> (1-3)])],
enable_debugging=$enableval, enable_debugging=no)

case $enable_debugging in
yes)  AC_DEFINE(kslDEBUGLEVEL, 0,[No debugging. ]);;
1)  AC_DEFINE(kslDEBUGLEVEL, 1,[Defines debugging level 1.]);;
2)  AC_DEFINE(kslDEBUGLEVEL, 2,[Defines debugging level 2.]);;
3)  AC_DEFINE(kslDEBUGLEVEL, 3,[Defines debugging level 3.]);;
no)  AC_DEFINE(kslDEBUGLEVEL, 0,[No debugging.]);;
*)  AC_MSG_ERROR([Unknown argument to --enable-debugging: $enable_debugging]);;
esac



if test "$enable_debugging" != "no"; then
AC_DEFINE(DEBUG,1,[Defines debugging .])
CFLAGS="$CFLAGS -ggdb  -Wall -m64 -std=gnu99    "
else
CFLAGS="$CFLAGS  -O2 -funroll-loops   -Wall  -std=gnu99 "
DEBUG=0
fi


AC_ARG_ENABLE(valgrind-tests,
[--enable-valgrind-tests run self tests under valgrind],
opt_valgrind_tests="$enableval", opt_valgrind_tests="no")

AC_ARG_ENABLE(trace,
[--enable-trace run trace tests under valgrind],
opt_prof_tests="$enableval", opt_prof_tests="no")

TRACEOPTS="--tool=callgrind "

OPTS="-q --trace-children=yes --error-exitcode=1 --dsymutil=yes --leak-check=yes --db-attach=yes "

VALGRIND=

# Run self-tests under valgrind?
if test "$opt_valgrind_tests" = "yes"; then
AC_DEFINE(DEBUG,1,[Defines debugging .])
CFLAGS=" -ggdb -pedantic -Wall -m64 -std=gnu99"
AC_CHECK_PROGS(VALGRIND, valgrind)
VALGRIND="valgrind $OPTS"
fi

if test "$opt_prof_tests" = "yes"; then
CFLAGS="-DDEBUG -ggdb -pedantic -Wall -m64 -std=gnu99"
AC_CHECK_PROGS(VALGRIND, valgrind)
VALGRIND="valgrind $TRACEOPTS"
fi


AC_MSG_CHECKING([whether self tests are run under valgrind])
AC_MSG_RESULT($opt_valgrind_tests)
AC_MSG_RESULT($opt_prof_tests)


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
AC_PROG_RANLIB

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
AC_CONFIG_HEADER([src/config.h])
AC_CONFIG_FILES([lib/Makefile])

AC_OUTPUT


echo "
     Configuration:
     
     compiler: ${CFLAGS}
     valgrind: ${VALGRIND}

"



