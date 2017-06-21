.PHONY: all clean

SHELL=/usr/bin/env bash -euf -o pipefail

.SECONDARY:

all:
	snakemake

clean:


