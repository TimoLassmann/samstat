# Samstat 
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

SAMStat reports ength distribution, base quality distribution, mapping
statistics, mismatch, insertion and deletion error profiles. The output is a
single html page:

![Image of example output](https://user-images.githubusercontent.com/8110320/33819832-69b39698-de87-11e7-9dfd-bdc330bfe78c.jpg)


# Install

``` bash
git clone https://github.com/TimoLassmann/samstat.git
cd samstat 
./autogen.sh
make 
make check 
make install 
```

# Usage

``` sh
samstat <file.sam>  <file.bam>  <file.fa>  <file.fq> .... 
```

For each input file SAMStat will create a single html page named after the input file name plus a dot html suffix.

# Please cite:

Lassmann et al. (2010) "SAMStat: monitoring biases in next generation sequencing data." Bioinformatics doi:10.1093/bioinformatics/btq614 [PMID: 21088025] 
[doi:10.1093/bioinformatics/btq614](http://dx.doi.org/10.1093%2Fbioinformatics%2Fbtq614)

