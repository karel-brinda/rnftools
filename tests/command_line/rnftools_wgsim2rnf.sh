#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

eval FA="humanmito.fa"
eval WGSIM="~/.smbl/bin/wgsim"

$WGSIM \
	-N 100 \
	-1 100 \
	-2 100 \
	$FA wgsim_1.fq wgsim_2.fq \
	> /dev/null


# 1) SE test, no contamination
rnftools wgsim2rnf \
	--faidx ${FA}.fai \
	--wgsim-fastq-1 wgsim_1.fq \
	--rnf-fastq _wgsim_rnf_se.fq \


# 2) PE test, no contamination
rnftools wgsim2rnf \
	--faidx ${FA}.fai \
	--wgsim-fastq-1 wgsim_1.fq \
	--wgsim-fastq-2 wgsim_2.fq \
	--rnf-fastq _wgsim_rnf_pe.fq \
