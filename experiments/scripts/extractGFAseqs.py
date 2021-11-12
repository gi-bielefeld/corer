#!/usr/bin/env python3

from sys import argv

counter = 1
K = 21

for l in open(argv[1], 'r'):
	if l.startswith('S'):
		print('>%d' %counter)
		print(l.split('\t')[2])
		counter += 1
		