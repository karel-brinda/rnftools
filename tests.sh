#! /usr/bin/env bash

set -e -o pipefail; 

(
	echo "Example 1"
	cd examples/01_simple_read_simulation
	snakemake
)

(
	echo "Example 2"
	cd examples/02_more_complex_read_simulation
	snakemake
)

(
	echo "Example 3"
	cd examples/03_mapper_evaluation_SE
	
	cd bams
	snakemake
	cd ..

	snakemake
)




