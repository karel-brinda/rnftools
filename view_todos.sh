#! /usr/bin/env bash

find rnftools -name "*.py" -type f | xargs grep -ni --color=auto TODO

