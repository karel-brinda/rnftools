#! /usr/bin/env bash

set -e

cd "$(dirname "$0")"

(
	rm -fr bats
	git clone git://github.com/sstephenson/bats
	cd bats
	./install.sh /usr/local
)

rm -fr bats
