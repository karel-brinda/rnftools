#! /usr/bin/env python3

exec(open("rnftools/version.py").read())

numbers=VERSION.split(".")
numbers[-1]=str(int(numbers[-1])+1)

version=".".join(numbers)

with open("rnftools/version.py","w") as f:
	f.write('VERSION="{}"'.format(version))
