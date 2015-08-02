#! /usr/bin/env bash

set -o pipefail
#set -x

(
	cd $1

	pwd
	for f in *.fq;
	do
		(
			echo "$f"
			rnftools validate -w -i $f
		)
	done;
)
