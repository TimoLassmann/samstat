
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@ 

MYLIBDIRS = ../tldevel/libthrpool.a \
../tldevel/libtldevel.a 


MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
-I${top_srcdir}/${LIB_TLDEVELDIR}

AM_CPPFLAGS=${MYINCDIRS} 

bin_PROGRAMS = samstat



samstat_SOURCES = \
samstat.h \
interface.c \
nuc_code.h \
nuc_code.c \
misc.h \
misc.c \
main.c \
io.c \
io.h \
hmm.c \
hmm.h \
viz.h \
viz.c


samstat_LDADD = ${MYLIBDIRS}

TESTS = hmm_ITEST

TESTS_ENVIRONMENT = \
	devdir='$(top_srcdir)'/dev; \
         export devdir; \
	$(VALGRIND)

check_PROGRAMS = hmm_ITEST


hmm_ITEST_SOURCES = hmm.c hmm.h
hmm_ITEST_LDADD = ${MYLIBDIRS}
hmm_ITEST_CPPFLAGS  = ${MYINCDIRS}  -DITEST



