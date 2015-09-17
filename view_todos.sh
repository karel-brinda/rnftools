#! /usr/bin/env bash

find rnftools -name "*.py" -type f | xargs grep TODO
find rnftools -name "*.py" -type f | xargs grep FIXME
