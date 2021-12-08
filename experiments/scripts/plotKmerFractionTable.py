#!/usr/bin/env python3

import sys

#This script plots an overview table of all prokaryotic datasets and their core k-mer fractions

SPECIES = ["bifidobacteriumAnimalis", "yersiniaPestis", "enterococcusFaecium", "listeriaMonocytogenes"]
PANGENOME_SIZES = [18, 48, 153, 263]

fractions = {}

for f in sys.argv[1:]:
	for s in SPECIES:
		elems = open(f, 'r').readline().split(' ')

		if f.find(s) >= 0:
			fractions[s] = float(elems[8]) / float(elems[10])

print('\t'.join(["species", "n", r"core $k$-mer fraction"]))

for i in range(len(SPECIES)):
	print('\t'.join([SPECIES[i], str(PANGENOME_SIZES[i]), "%.2f" %fractions[SPECIES[i]]]))
	