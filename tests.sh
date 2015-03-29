#! /usr/bin/env bash

set -e

echo
echo
echo
echo
echo
echo
echo "===================================================================="
echo
echo
echo
echo "                 TEST: example 1"
echo
echo
echo
echo "===================================================================="
echo
echo
echo

(
	set -ex -o pipefail; 

	cd examples/01_simple_read_simulation
	snakemake -p --cores
)

echo
echo
echo
echo
echo
echo
echo "===================================================================="
echo
echo
echo
echo "                 TEST: example 2" 
echo
echo
echo
echo "===================================================================="
echo
echo
echo

(
	set -ex -o pipefail; 

	cd examples/02_more_complex_read_simulation
	snakemake -p --cores
)

echo
echo
echo
echo
echo
echo
echo "===================================================================="
echo
echo
echo
echo "                 TEST: example 3" 
echo
echo
echo
echo "===================================================================="
echo
echo
echo

(
	set -ex -o pipefail; 

	cd examples/03_mapper_evaluation_SE
	
	cd bams
	snakemake -p --cores
	cd ..

	cd report
	snakemake -p --cores
)




echo
echo
echo
echo
echo
echo
echo "===================================================================="
echo
echo
echo
echo "                 TEST: example 4" 
echo
echo
echo
echo "===================================================================="
echo
echo
echo

(
	set -ex -o pipefail; 

	cd examples/04_mapper_evaluation_PE
	
	cd bams
	snakemake -p --cores
	cd ..

	cd report
	snakemake -p --cores
)




