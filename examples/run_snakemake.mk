.PHONY: all clean

SHELL=/usr/bin/env bash
.SHELLFLAGS = -euf -o pipefail

.SECONDARY:

all:
	snakemake

clean:


