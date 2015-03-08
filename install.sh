#! /usr/bin/env bash

set -eu -o pipefail

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

