#! /usr/bin/env bash

set -eu
set -o verbose
set -o pipefail

cd "$(dirname "$0")"

rnftools check
