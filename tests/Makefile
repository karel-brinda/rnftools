.PHONY: all clean

SHELL=/usr/bin/env bash
.SHELLFLAGS=-e -c -o pipefail

SUBDIRS = $(wildcard */)

default: all

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

all clean : $(SUBDIRS)

