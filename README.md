




Next generation sequencing is being applied to understand individual variation,
the RNA output of a cell and epigenetic regulation. The millions of sequenced
reads are commonly stored in fasta, fastq and after mapping to a reference
genome in the alignment / map format (SAM/BAM). To monitor the sequence quality
over time and to identify problems it is necessary to report various statistics
of the reads at different stages during processing.

SAMStat is an efficient C program to quickly display statistics of large
sequence files from next generation sequencing projects. When applied to SAM/BAM
files all statistics are reported for unmapped, poorly and accurately mapped
reads separately. This allows for identification of a variety of problems, such
as remaining linker and adaptor sequences, causing poor mapping. Apart from this
SAMStat can be used to verify individual processing steps in large analysis
pipelines.

SAMStat reports nucleotide composition, length distribution, base quality
distribution, mapping statistics, mismatch, insertion and deletion error
profiles, di-nucleotide and 10-mer over-representation. The output is a single
html5 page which can be interpreted by a non-specialist.


README.md

SAMstat version 1.5

Installation:

Unpack the tarball:
bash-3.1$ tar -zxvf samstat-XXX.tar.gz
bash-3.1$ cd samstat
bash-3.1$ ./configure

bash-3.1$ make
bash-3.1$ make check

At this point the samstat executable appears in the src directory. You can copy it to any directory in your path. To install it system wide type:

bash-3.1$ make install
