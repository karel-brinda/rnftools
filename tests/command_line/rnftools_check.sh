#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

rnftools check
