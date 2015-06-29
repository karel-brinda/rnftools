#! /usr/bin/env bash

for x in check publication sam2rnf art2rnf curesim2rnf dwgsim2rnf mason2rnf wgsim2rnf merge sam2es es2et et2roc sam2roc ; do
	fn="$x.txt"
	echo "$fn"
	echo "$ rnftools $x -h" > $fn
	echo >> $fn
	rnftools $x -h >> $fn
done
