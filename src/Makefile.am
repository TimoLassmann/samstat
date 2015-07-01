
AM_CPPFLAGS = -I$(top_srcdir)/lib

bin_PROGRAMS = samstat



samstat_SOURCES = \
malloc_macro.h \
interface.h \
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

samstat_LDADD = -lm -lpthread  $(abs_top_builddir)/lib/libks.a

TESTS = hmm_ITEST

TESTS_ENVIRONMENT = \
	devdir='$(top_srcdir)'/dev; \
         export devdir; \
	$(VALGRIND)

check_PROGRAMS = hmm_ITEST


hmm_ITEST_SOURCES = hmm.c hmm.h
hmm_ITEST_LDADD = $(samstat_LDADD)
hmm_ITEST_CPPFLAGS  =  -DITEST



