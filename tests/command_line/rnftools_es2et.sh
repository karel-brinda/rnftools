#! /usr/bin/env bash

set -eu
set -o verbose
set -o pipefail

cd "$(dirname "$0")"

rnftools es2et \
	-i rnf_alignment.es \
	-o - > /dev/null
