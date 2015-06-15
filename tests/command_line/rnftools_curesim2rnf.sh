#! /usr/bin/env bash

set -eux
set -o pipefail

eval FA="~/.smbl/fa/Mycobacterium_tuberculosis.fa"
eval CURESIM_JAR="~/.smbl/bin/curesim.jar"
CURESIM="java -jar $CURESIM_JAR"

$CURESIM \
	-n 100 \
	-f $FA \
	-o curesim.fq \

rnftools curesim2rnf \
	--fasta-index ${FA}.fai \
	--curesim-fastq curesim.fq \
	--rnf-fastq _curesim_rnf_se.fq \


