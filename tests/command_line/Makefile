SCRIPTS=$(wildcard *.sh)

.PHONY: all clean $(SCRIPTS)

.SECONDARY:

SHELL=/usr/bin/env bash -e -c -o pipefail


all: $(SCRIPTS)

$(SCRIPTS) :
	echo "script" $<
	bash $@

clean:


