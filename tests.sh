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
echo "                 TESTS: command line"
echo
echo
echo
echo "===================================================================="
echo
echo
echo

(
	set -e -o pipefail; 

	cd tests/command_line
	snakemake -p -s __ensure_programs.snake --cores
	for e in *.sh; do
		(
			echo "--------------------------------------------------------------------"
			echo
			echo "                 TEST: $e"
			echo
			echo "--------------------------------------------------------------------"

			bash $e
		)
	done
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
echo "                 TESTS: basic snakemake tests"
echo
echo
echo
echo "===================================================================="
echo
echo
echo

(
	set -e -o pipefail; 

	for e in tests/snakemake/*; do
		(
			cd "$e"

			echo "--------------------------------------------------------------------"
			echo
			echo "                 TEST: $e"
			echo
			echo "--------------------------------------------------------------------"


			snakemake -p --cores
		)
	done
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
echo "                 TESTS: tutorial examples"
echo
echo
echo
echo "===================================================================="
echo
echo
echo

(
	set -e -o pipefail; 

	for e in examples/01_tutorial/02_simulation; do
		(
			cd "$e"

			for d in *; do
				(
					echo "--------------------------------------------------------------------"
					echo
					echo "                 TEST: $d"
					echo
					echo "--------------------------------------------------------------------"


					cd "$d"
					snakemake -p --cores
				)
			done
		)
	done

	for e in examples/01_tutorial/03_evaluation; do
		(
			cd "$e"

			for d in *; do
				(
					echo "--------------------------------------------------------------------"
					echo
					echo "                 TEST: $d"
					echo
					echo "--------------------------------------------------------------------"


					cd "$d"
					(
						cd bams && snakemake -p --cores
					)
					(
						cd report && snakemake -p --cores
					)
				)
			done
		)
	done
)
