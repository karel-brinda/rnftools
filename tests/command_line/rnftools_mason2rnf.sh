#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

eval FA="humanmito.fa"
eval MASON="~/.smbl/bin/mason_simulator"

# 1) SE test, no contamination
$MASON \
	--input-reference $FA \
	--illumina-read-length 100 \
	--num-fragments 100 \
	--out-alignment mason_se.sam \
	--out tmp.mason.1.fq \

rnftools mason2rnf \
	--faidx ${FA}.fai \
	--sam mason_se.sam \
	--rnf-fastq _mason_rnf_se.fq \

# 2) PE test, no contamination
$MASON \
	--input-reference $FA \
	--illumina-read-length 100 \
	--num-fragments 100 \
	--out-alignment mason_pe.sam \
	--out tmp.mason.1.fq \
	--out-right tmp.mason.2.fq \

rnftools mason2rnf \
	--faidx ${FA}.fai \
	--sam mason_pe.sam \
	--rnf-fastq _mason_rnf_pe.fq \
