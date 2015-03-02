#! /usr/bin/env bash

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
	rm -fR ~/.smbl

	cd examples/01_simple_read_simulation
	snakemake -p
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
	rm -fR ~/.smbl

	cd examples/02_more_complex_read_simulation
	snakemake -p
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
	rm -fR ~/.smbl

	cd examples/03_mapper_evaluation_SE
	
	cd bams
	snakemake -p
	cd ..

	snakemake -p
)




