.PHONY: all clean

SHELL=/usr/bin/env bash
.SHELLFLAGS = -euf -o pipefail

.SECONDARY:

$(info )
$(info )
$(info )
$(info )
$(info )
$(info =======================================================)
$(info $(shell pwd))
$(info =======================================================)
$(info )

all:
	snakemake

clean:


