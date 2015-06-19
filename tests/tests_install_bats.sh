#! /usr/bin/env bash

set -e

cd "$(dirname "$0")"

(
	git clone git://github.com/sstephenson/bats
	cd bats
	./install.sh /usr/local
)

rm -fr bats
