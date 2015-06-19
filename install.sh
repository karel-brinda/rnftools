#! /usr/bin/env bash

set -eu -o pipefail

cd "$(dirname "$0")"

rm -fR build dist RNFtools.egg-info

echo
echo
echo "   TEST"
echo
echo

python3 setup.py check

echo
echo
echo "   INSTALLATION"
echo
echo

python3 setup.py install

