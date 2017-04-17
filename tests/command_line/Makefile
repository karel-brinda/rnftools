SCRIPTS=$(wildcard *.sh)

.PHONY: all clean $(SCRIPTS)

.SECONDARY:

SHELL=/usr/bin/env bash
.SHELLFLAGS = -e -c -o pipefail


all: $(SCRIPTS)

$(SCRIPTS) :
	echo "script" $<
	bash $@

clean:


