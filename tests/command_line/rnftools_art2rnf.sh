#! /usr/bin/env bash

set -eux
set -o pipefail

eval FA="humanmito.fa"
eval ART_ILLUMINA="~/.smbl/bin/art_illumina"


# 1) SE test, no contamination
$ART_ILLUMINA -sam \
	--in $FA \
	--len 100 \
	--fcov 1 \
	--out art

rnftools art2rnf \
	--fasta-index ${FA}.fai \
	--sam art.sam \
	--fastq _art_rnf_se.fq \

# 2) PE test, no contamination
#rnftools art2rnf --fasta-index ${FA}.fai --fastq reads_in_rnf.fq --wgsim-fastq-1 wgsim_1.fq --wgsim-fastq-2 wgsim_2.fq
