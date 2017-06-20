.PHONY: all clean

SHELL=/usr/bin/env bash -eufc -o pipefail

.SECONDARY:

all: bwamem.bam

bwamem.bam: ../1_reads/SE_reads.fq
	bwa index ../reference.fa
	bwa mem ../reference.fa ../1_reads/SE_reads.fq | samtools view -b > $@

clean:
	rm -f *.bam
