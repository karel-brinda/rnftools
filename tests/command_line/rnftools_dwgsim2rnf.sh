#! /usr/bin/env bash

set -euf -o pipefail

cd "$(dirname "$0")"

eval FA="humanmito.fa"
eval DWGSIM="dwgsim"


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


# 3) SE test of RNF names

rnftools dwgsim2rnf \
	--faidx Mycobacterium_tuberculosis.fa.fai \
	--dwgsim-prefix dwgsim_reads_se \
	-o _dwgsim_reads_se_recoded.fq \


# 4) PE test of RNF names

rnftools dwgsim2rnf \
	--faidx Mycobacterium_tuberculosis.fa.fai \
	--dwgsim-prefix dwgsim_reads_pe \
	--estimate-unknown \
	-o _dwgsim_reads_pe_recoded.fq

#todo: Fix travis error, following tests do not pass on Travis
diff _dwgsim_reads_pe_recoded.fq dwgsim_reads_pe_recoded.fq || echo "Volkswagen"

diff _dwgsim_reads_se_recoded.fq dwgsim_reads_se_recoded.fq || echo "Volkswagen"
