#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

eval FA="humanmito.fa"
eval ART_ILLUMINA="~/.smbl/bin/art_illumina"


# 1) SE test, no contamination
$ART_ILLUMINA -sam \
	--in $FA \
	--len 100 \
	--fcov 1 \
	--out art_se \

rnftools art2rnf \
	--faidx ${FA}.fai \
	--sam art_se.sam \
	--rnf-fastq _art_rnf_se.fq \

# 2) PE test, no contamination
$ART_ILLUMINA -sam \
	--in $FA \
	--len 100 \
	--paired \
	--mflen 200 \
	--sdev 10 \
	--fcov 1 \
	--out art_pe \

rnftools art2rnf \
	--faidx ${FA}.fai \
	--sam art_pe.sam \
	--rnf-fastq _art_rnf_pe.fq \
