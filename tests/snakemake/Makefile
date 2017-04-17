DIRS=$(wildcard */run.sh)

.PHONY: all clean $(DIRS)

.SECONDARY:

SHELL=/usr/bin/env bash
.SHELLFLAGS = -e -c -o pipefail


all: $(DIRS)

$(DIRS) :
	bash $@

clean:


