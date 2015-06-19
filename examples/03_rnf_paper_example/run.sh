#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

(
	cd bams
	snakemake "$@"
)

(
	cd reports
	snakemake "$@"
)
