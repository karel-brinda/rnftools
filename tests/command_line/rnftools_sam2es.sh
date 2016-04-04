#! /usr/bin/env bash

set -eu
set -o verbose
set -o pipefail

cd "$(dirname "$0")"

# SAM
rnftools sam2es \
	-i rnf_alignment.sam \
	-o - > /dev/null

# BAM
rnftools sam2es \
	-i rnf_alignment.bam \
	-o - > /dev/null