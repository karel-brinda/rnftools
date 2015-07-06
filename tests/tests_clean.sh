#! /usr/bin/env bash

set -e

for d in . ../examples;
do
	(
		cd $d
		git clean -fdx
	)
done;

