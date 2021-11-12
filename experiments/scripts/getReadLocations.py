#!/usr/bin/env python3

import sys

isFirst = True
addresses = []

for l in open(sys.argv[1]):
	if isFirst:
		isFirst = False
		continue

	for a in l.split('\t')[6].split(';'):
		addresses.append(a)

print(' '.join(addresses))