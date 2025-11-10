.PHONY: all clean

SHELL=/usr/bin/env bash -euf -o pipefail
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
	snakemake --cores all -p

clean:


