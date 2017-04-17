#! /usr/bin/env bash

set -euf -o pipefail

cd "$(dirname "$0")"

rnftools check
