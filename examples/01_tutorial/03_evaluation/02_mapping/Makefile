.PHONY: all clean

SHELL=/usr/bin/env bash -eufc -o pipefail

.SECONDARY:

all: bwa-mem.bam bwa-sw.bam

bwa-mem.bam: ../01_reads/SE_reads.fq
	bwa index ../reference.fa
	bwa mem ../reference.fa ../01_reads/SE_reads.fq | samtools view -b > $@

bwa-sw.bam: ../01_reads/SE_reads.fq
	bwa index ../reference.fa
	bwa bwasw ../reference.fa ../01_reads/SE_reads.fq | samtools view -b > $@

clean:
	rm -f *.bam
