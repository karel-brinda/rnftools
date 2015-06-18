#! /usr/bin/env bash

set -eux
set -o pipefail

rnftools et2roc \
	-i rnf_alignment.et \
	-o - > /dev/null
