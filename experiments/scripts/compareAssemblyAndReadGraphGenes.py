#!/usr/bin/env python3

import argparse as args

#This script compares two cores predicted by Corer on the gene level

def parseBifrostRes(resFiles):
	coreGenes = {}

	for f in resFiles:
		isFrst = True
		accName = f.split("_q")[1].split(".tsv")[0]

		for l in open(f, 'r'):
			if isFrst:
				isFrst = False
				continue

			gene, flag = l.split('\t')

			if int(flag) > 0:
				coreGenes[accName + gene] = None

	return coreGenes

if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script compares two cores predicted by Corer on the gene level.")
	parser.add_argument('-c1', metavar='CoreGenes1', type=str, required=True, nargs='+', help="Bifrost query results for first core")
	parser.add_argument('-c2', metavar='CoreGenes2', type=str, required=True, nargs='+', help="Bifrost query results for second core")

	arguments = parser.parse_args()

	#Parse results
	corerGenes1 = parseBifrostRes(arguments.c1)
	corerGenes2 = parseBifrostRes(arguments.c2)

	firstCoreOnly = 0
	secondCoreOnly = 0
	shared = 0

	for g in corerGenes1:
		if g in corerGenes2:
			shared += 1
		else:
			firstCoreOnly += 1

	for g in corerGenes2:
		if not g in corerGenes1:
			secondCoreOnly += 1

	print("%d genes are shared" %shared)
	print("%d genes are exclusively present in first core" %firstCoreOnly)
	print("%d genes are exclusively present in second core" %secondCoreOnly)
