#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

eval FA="humanmito.fa"
eval DWGSIM="~/.smbl/bin/dwgsim"


# 1) SE test, no contamination

$DWGSIM \
	-y 0 \
	-N 100 \
	-1 100 \
	-2 000 \
	$FA sim_se

rnftools dwgsim2rnf \
	--faidx ${FA}.fai \
	--dwgsim-prefix sim_se \
	--rnf-fastq _dwgsim_rnf_se.fq \


# 2) PE test, no contamination

$DWGSIM \
	-y 0 \
	-N 100 \
	-1 100 \
	-2 100 \
	$FA sim_pe \

rnftools dwgsim2rnf \
	--faidx ${FA}.fai \
	--dwgsim-prefix sim_pe \
	--rnf-fastq _dwgsim_rnf_pe.fq \
