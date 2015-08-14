#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

eval FA="humanmito.fa"
eval CURESIM_JAR="~/.smbl/bin/CuReSim.jar"
CURESIM="java -jar $CURESIM_JAR"

$CURESIM \
	-n 100 \
	-f $FA \
	-o curesim.fq \
2> /dev/null

rnftools curesim2rnf \
	--faidx ${FA}.fai \
	--curesim-fastq curesim.fq \
	--rnf-fastq _curesim_rnf_se.fq \
