#! /usr/bin/env bash

set -eux
set -o pipefail

eval FA="~/.smbl/fa/Mycobacterium_tuberculosis.fa"
eval DWGSIM="~/.smbl/bin/dwgsim"


# 1) SE test, no contamination

$DWGSIM \
	-y 0 \
	-N 100 \
	-1 100 \
	-2 000 \
	$FA sim_se

rnftools dwgsim2rnf \
	--fasta-index ${FA}.fai \
	--dwgsim-prefix sim_se \
	--fastq _dwgsim_rnf_se.fq \


# 2) PE test, no contamination

$DWGSIM \
	-y 0 \
	-N 100 \
	-1 100 \
	-2 100 \
	$FA sim_pe \

rnftools dwgsim2rnf \
	--fasta-index ${FA}.fai \
	--dwgsim-prefix sim_pe \
	--fastq _dwgsim_rnf_pe.fq \
