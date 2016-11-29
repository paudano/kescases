ALL_TARGETS=bin/time bin/traceproc bin/bwa bin/samtools bin/spades.py bin/python bin/snakemake bin/tabix bin/rtg

.PHONY: all
all: $(ALL_TARGETS)

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

