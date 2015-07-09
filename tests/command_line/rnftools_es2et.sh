#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

rnftools es2et \
	-i rnf_alignment.es \
	-o - > /dev/null
