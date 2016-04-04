#! /usr/bin/env bash

set -eu
set -o verbose
set -o pipefail

cd "$(dirname "$0")"

rnftools et2roc \
	-i rnf_alignment.et \
	-o - > /dev/null
