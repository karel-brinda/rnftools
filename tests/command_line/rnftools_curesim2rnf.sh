#! /usr/bin/env bash

set -euf -o pipefail

cd "$(dirname "$0")"

eval FA="humanmito.fa"
CURESIM="curesim"

$CURESIM \
	-n 100 \
	-f $FA \
	-o curesim.fq \
2> /dev/null

rnftools curesim2rnf \
	--faidx ${FA}.fai \
	--curesim-fastq curesim.fq \
	--rnf-fastq _curesim_rnf_se.fq \
