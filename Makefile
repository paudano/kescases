ALL_FREE_TARGETS=bin/time bin/traceproc bin/bwa bin/samtools bin/spades.py bin/python bin/snakemake bin/tabix bin/rtg bin/R bin/Rscript bin/fastq-dump bin/blastn bin/quast.py

ALL_NONFREE_TARGETS=lib/GenomeAnalysisTK.jar lib/picard.jar

.PHONY: all
all: $(ALL_NONFREE_TARGETS) $(ALL_FREE_TARGETS)

.PHONY: allfree
allfree: $(ALL_FREE_TARGETS)

lib/GenomeAnalysisTK.jar lib/picard.jar:
	@ if [ ! -a $@ ]; then \
		echo "* $@ is not installed and must be manually downloaded." ; \
		echo "* See lib/NOTES for installation instructions." ; \
		echo "* To build all freely-available software, run \"make allfree\"" ; \
		exit 1 ; \
	fi ; \
	echo OK

bin/quast.py:
	make -C build/quast

bin/blastn:
	make -C build/blast

bin/fastq-dump:
	make -C build/sra

bin/R bin/Rscript:
	make -C build/r

bin/rtg:
	make -C build/rtg

bin/tabix:
	make -C build/tabix

bin/python:
	make -C build/miniconda ../../bin/python

bin/snakemake:
	make -C build/miniconda ../../bin/snakemake

bin/time:
	make -C build/time

bin/traceproc:
	make -C build/traceproc

bin/bwa:
	make -C build/bwa

bin/samtools:
	make -C build/samtools

bin/spades.py:
	make -C build/spades

.PHONY: clean
clean:
	make -C build/time clean
	make -C build/traceproc clean
	make -C build/bwa clean
	make -C build/samtools clean
	make -C build/spades clean

