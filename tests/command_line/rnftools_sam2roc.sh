#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

# SAM
rnftools sam2roc \
	-i rnf_alignment.sam \
	-o - > /dev/null

# BAM
rnftools sam2roc \
	-i rnf_alignment.bam \
	-o - > /dev/null