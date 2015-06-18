#! /usr/bin/env bash

set -eux
set -o pipefail

rnftools es2et \
	-i rnf_alignment.es \
	-o - > /dev/null
