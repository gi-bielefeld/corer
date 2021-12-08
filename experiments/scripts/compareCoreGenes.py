#!/usr/bin/env python3

import argparse as args

#This script compares the cores of Corer and SibeliaZ on the gene level

#Setting up the argument parser
parser = args.ArgumentParser(description="This script compares the cores of Corer and SibeliaZ on the gene level.")
parser.add_argument('-c', metavar='CorerGenes', type=str, required=True, nargs='+', help="Bifrost query results")
parser.add_argument('-s', metavar='SibeliaZGenes', type=str, required=True, nargs='+', help="Gene coverage files")
parser.add_argument('-t', metavar='covThres', type=float, required=True, help="Gene coverage threshold")

arguments = parser.parse_args()

#Parse Corer's results
corerGenes = {}

for f in arguments.c:
	isFrst = True
	accName = f.split("_q")[1].split(".tsv")[0]

	for l in open(f, 'r'):
		if isFrst:
			isFrst = False
			continue

		gene, flag = l.split('\t')

		if int(flag) > 0:
			corerGenes[accName + gene] = None

#Parse SibeliaZ's results
sibeliaGenes = {}
	    
for f in arguments.s:
	accName = f.split("_q")[1].split(".txt")[0]

	for l in open(f, 'r'):
		if len(l.split(' ')) < 3:
			l = l.strip()
			gName, relCov = l.split(' ')

			if float(relCov) >= arguments.t:
				sibeliaGenes[accName + gName] = float(relCov)

corerOnly = 0
sibeliaOnly = 0
shared = 0

for g in corerGenes:
	if g in sibeliaGenes:
		shared += 1
	else:
		corerOnly += 1

for g in sibeliaGenes:
	if not g in corerGenes:
		sibeliaOnly += 1

print("%d genes are shared" %shared)
print("%d genes are exclusively found by Corer" %corerOnly)
print("%d genes are exclusively found by SibeliaZ" %sibeliaOnly)
