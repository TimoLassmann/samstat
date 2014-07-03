
VFLAGS = --track-origins=yes --tool=memcheck --leak-check=yes --error-exitcode=1


AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}



SUBDIRS = src 

dist_noinst_SCRIPTS = autogen.sh




rsync:
	ssh lm3 'rm -rf lm3:/data/lassmann/tmp/testing '
	rsync -ra -e ssh --exclude '/.git' --exclude '/bin' --exclude '/doc' --exclude '*.tar.gz' --exclude '/build' . lm3:/data/lassmann/tmp/testing

val: rsync
	ssh lm3 'cd /data/lassmann/tmp/testing && make clean && make check '


benchmark:
	


