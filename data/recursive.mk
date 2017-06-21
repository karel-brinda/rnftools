.PHONY: all clean

SHELL=/usr/bin/env bash -e -c -o pipefail

SUBDIRS = $(wildcard */)

default: all

$(info )
$(info )
$(info )
$(info )
$(info )
$(info =======================================================)
$(info $(shell pwd))
$(info =======================================================)
$(info )

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

all clean : $(SUBDIRS)

