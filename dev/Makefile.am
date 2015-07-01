AUTOMAKE_OPTIONS = serial-tests

HASHMARK = \#

EXTRA_DIST= run_samstat_tests.sh aux.fa c1.fa ce.fa xx.fa aux$(HASHMARK)aux.sam c1$(HASHMARK)pad3.sam ce$(HASHMARK)large_seq.sam ce$(HASHMARK)unmap2.sam xx$(HASHMARK)minimal.sam c1$(HASHMARK)bounds.sam ce$(HASHMARK)1.sam ce$(HASHMARK)tag_depadded.sam fieldarith.sam xx$(HASHMARK)pair.sam c1$(HASHMARK)clip.sam ce$(HASHMARK)2.sam ce$(HASHMARK)tag_padded.sam xx$(HASHMARK)blank.sam xx$(HASHMARK)rg.sam c1$(HASHMARK)pad1.sam ce$(HASHMARK)5.sam ce$(HASHMARK)unmap.sam xx$(HASHMARK)large_aux.sam xx$(HASHMARK)triplet.sam c1$(HASHMARK)pad2.sam ce$(HASHMARK)5b.sam ce$(HASHMARK)unmap1.sam xx$(HASHMARK)large_aux2.sam xx$(HASHMARK)unsorted.sam



TESTS_ENVIRONMENT = \
testdatafiledir='$(top_srcdir)'/dev; \
export testdatafiledir;
	$(VALGRIND)

TESTS = run_samstat_tests.sh




CLEANFILES = *html



DISTCLEANFILES = *html