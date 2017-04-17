#! /usr/bin/env bash

set -euf -o pipefail

cd "$(dirname "$0")"

rnftools es2et \
	-i rnf_alignment.es \
	-o - > /dev/null
