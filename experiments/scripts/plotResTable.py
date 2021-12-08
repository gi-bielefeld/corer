#!/usr/bin/env python3

import argparse as args
from math import ceil
from copy import deepcopy

#This script takes the results of Corer, Panaroo and SibeliaZ, parses them with respect to a given core quorum and gene coverage threshold, and outputs the results in tabular format.

def calcIntSec(aRes, bRes):
	return len([g for g in aRes if g in bRes])

def calcDiff(aRes, bRes):
	return len([g for g in aRes if not g in bRes])

def calcPerc(common, exclusive):
	return  (common / float(common + exclusive)) * 100

if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script outputs the given results of Corer, Panaroo and SibeliaZ in tabular format.")
	parser.add_argument('-b', metavar='BifrostRes', type=str, required=True, nargs='+', help="Bifrost query results")
	parser.add_argument('-p', metavar='PanarooRes', type=str, required=True, nargs='+', help="Panaroo results")
	parser.add_argument('-s', metavar='SibeliaZRes', type=str, required=True, nargs='+', help="SibeliaZ results")
	parser.add_argument('-q', metavar='relQuorum', type=float, required=True, help="Relative core quorum")
	parser.add_argument('-t', metavar='covThres', type=float, required=True, help="Gene coverage threshold")

	arguments = parser.parse_args()
	toolsRes = []

	#Parse Corer's results
	toolsRes.append({})

	for f in arguments.b:
		isFrst = True
		species = f.split('/')[4]
		accName = f.split("_q")[1].split(".tsv")[0]

		if not species in toolsRes[0]:
			toolsRes[0][species] = {}

		for l in open(f, 'r'):
			if isFrst:
				isFrst = False
				continue

			gene, flag = l.split('\t')

			if int(flag) > 0:
				toolsRes[0][species][accName + gene] = None

	#Parse Panaroo's results
	toolsRes.append({})

	for f in arguments.p:
		isFrst = True
		species = f.split('/')[4]

		if not species in toolsRes[1]:
			toolsRes[1][species] = {}

		for l in open(f, 'r'):
			l = l.strip()
			columns = l.split(',')[3:]
			panSize = len(columns)

			if isFrst:
				accNames = columns
				isFrst = False
				continue

			if columns.count('') <= panSize - ceil(arguments.q * panSize):
				for k in [accNames[i] + columns[i] for i in range(len(accNames)) if columns[i] != '' and columns[i].find("refound") < 0]:
					if not k in toolsRes[1][species]:
						toolsRes[1][species][k] = None

	#Parse SibeliaZ's results
	toolsRes.append({})
	    
	for f in arguments.s:
		species = f.split('/')[4]
		accName = f.split("_q")[1].split(".txt")[0]

		if not species in toolsRes[2]:
			toolsRes[2][species] = {}

		for l in open(f, 'r'):
			if len(l.split(' ')) < 3:
				l = l.strip()
				gName, relCov = l.split(' ')

				if float(relCov) >= arguments.t:
					toolsRes[2][species][accName + gName] = float(relCov)

	#Calculate intersections
	iSecs = {}
	diffs = {}

	for s in toolsRes[0]:
		iSecs[s] = [[] for t in toolsRes]
		diffs[s] = deepcopy(iSecs[s])

		for i in range(len(toolsRes)):
			for j in range(len(toolsRes)):
				if i < j and s in toolsRes[i] and s in toolsRes[j]:
					iSecs[s][i].append(calcIntSec(toolsRes[i][s], toolsRes[j][s]))

				if i != j and s in toolsRes[i] and s in toolsRes[j]:
					diffs[s][i].append(calcDiff(toolsRes[i][s], toolsRes[j][s]))

	#Output table
	print("\t".join(["", "", "B. animalis", "", "", "Y. pestis", "", "", "E. faecium", "", "", "L. monocytogenes"]))
	print("\t".join(["", "Corer", "%d (%.1f%%)" %(iSecs["bifidobacteriumAnimalis"][0][0], calcPerc(iSecs["bifidobacteriumAnimalis"][0][0], diffs["bifidobacteriumAnimalis"][0][0])), \
		"%d (%.1f%%)" %(iSecs["bifidobacteriumAnimalis"][0][1], calcPerc(iSecs["bifidobacteriumAnimalis"][0][1], diffs["bifidobacteriumAnimalis"][0][1])), "Corer", "%d (%.1f%%)" \
		%(iSecs["yersiniaPestis"][0][0], calcPerc(iSecs["yersiniaPestis"][0][0], diffs["yersiniaPestis"][0][0])), "%d (%.1f%%)" %(iSecs["yersiniaPestis"][0][1], calcPerc(iSecs["yersiniaPestis"][0][1], \
		diffs["yersiniaPestis"][0][1])), "Corer", "%d (%.1f%%)" %(iSecs["enterococcusFaecium"][0][0], calcPerc(iSecs["enterococcusFaecium"][0][0], diffs["enterococcusFaecium"][0][0])), "%d (%.1f%%)" \
		%(iSecs["enterococcusFaecium"][0][1], calcPerc(iSecs["enterococcusFaecium"][0][1], diffs["enterococcusFaecium"][0][1])), "Corer", "%d (%.1f%%)" %(iSecs["listeriaMonocytogenes"][0][0], calcPerc( \
		iSecs["listeriaMonocytogenes"][0][0], diffs["listeriaMonocytogenes"][0][0]))]))
	print("\t".join([r"$i\\cap j$", "%d (%.1f%%)" %(iSecs["bifidobacteriumAnimalis"][0][0], calcPerc(iSecs["bifidobacteriumAnimalis"][0][0], diffs["bifidobacteriumAnimalis"][1][0])), "Panaroo", \
		"%d (%.1f%%)" %(iSecs["bifidobacteriumAnimalis"][1][0], calcPerc(iSecs["bifidobacteriumAnimalis"][1][0], diffs["bifidobacteriumAnimalis"][1][1])), "%d (%.1f%%)" %(iSecs["yersiniaPestis"][0][0], \
		calcPerc(iSecs["yersiniaPestis"][0][0], diffs["yersiniaPestis"][1][0])), "Panaroo", "%d (%.1f%%)" %(iSecs["yersiniaPestis"][1][0], calcPerc(iSecs["yersiniaPestis"][1][0], \
		diffs["yersiniaPestis"][1][1])), "%d (%.1f%%)" %(iSecs["enterococcusFaecium"][0][0], calcPerc(iSecs["enterococcusFaecium"][0][0], diffs["enterococcusFaecium"][1][0])), "Panaroo", "%d (%.1f%%)" \
		%(iSecs["enterococcusFaecium"][1][0], calcPerc(iSecs["enterococcusFaecium"][1][0], diffs["enterococcusFaecium"][1][1])), "%d (%.1f%%)" %(iSecs["listeriaMonocytogenes"][0][0], calcPerc( \
		iSecs["listeriaMonocytogenes"][0][0], diffs["listeriaMonocytogenes"][1][0])), "Panaroo"]))
	print("\t".join(["", "%d (%.1f%%)" %(iSecs["bifidobacteriumAnimalis"][0][1], calcPerc(iSecs["bifidobacteriumAnimalis"][0][1], diffs["bifidobacteriumAnimalis"][2][0])), "%d (%.1f%%)" \
		%(iSecs["bifidobacteriumAnimalis"][1][0], calcPerc(iSecs["bifidobacteriumAnimalis"][1][0], diffs["bifidobacteriumAnimalis"][2][1])), "SibeliaZ", "%d (%.1f%%)" %(iSecs["yersiniaPestis"][0][1], \
		calcPerc(iSecs["yersiniaPestis"][0][1], diffs["yersiniaPestis"][2][0])), "%d (%.1f%%)" %(iSecs["yersiniaPestis"][1][0], calcPerc(iSecs["yersiniaPestis"][1][0], diffs["yersiniaPestis"][2][1])), \
		"SibeliaZ", "%d (%.1f%%)" %(iSecs["enterococcusFaecium"][0][1], calcPerc(iSecs["enterococcusFaecium"][0][1], diffs["enterococcusFaecium"][2][0])), "%d (%.1f%%)" \
		%(iSecs["enterococcusFaecium"][1][0], calcPerc(iSecs["enterococcusFaecium"][1][0], diffs["enterococcusFaecium"][2][1])), "SibeliaZ", "", ""]))
	print("\t".join(["", "Corer", "%d (%.1f%%)" %(diffs["bifidobacteriumAnimalis"][0][0], calcPerc(diffs["bifidobacteriumAnimalis"][0][0], iSecs["bifidobacteriumAnimalis"][0][0])), "%d (%.1f%%)" \
		%(diffs["bifidobacteriumAnimalis"][0][1], calcPerc(diffs["bifidobacteriumAnimalis"][0][1], iSecs["bifidobacteriumAnimalis"][0][1])), "Corer", "%d (%.1f%%)" %(diffs["yersiniaPestis"][0][0], \
		calcPerc(diffs["yersiniaPestis"][0][0], iSecs["yersiniaPestis"][0][0])), "%d (%.1f%%)" %(diffs["yersiniaPestis"][0][1], calcPerc(diffs["yersiniaPestis"][0][1], iSecs["yersiniaPestis"][0][1])), \
		"Corer", "%d (%.1f%%)" %(diffs["enterococcusFaecium"][0][0], calcPerc(diffs["enterococcusFaecium"][0][0], iSecs["enterococcusFaecium"][0][0])), "%d (%.1f%%)" %(diffs["enterococcusFaecium"][0][1], \
		calcPerc(diffs["enterococcusFaecium"][0][1], iSecs["enterococcusFaecium"][0][1])), "Corer", "%d (%.1f%%)" %(diffs["listeriaMonocytogenes"][0][0], calcPerc(diffs["listeriaMonocytogenes"][0][0], \
		iSecs["listeriaMonocytogenes"][0][0]))]))
	print("\t".join([r"$i\\backslash j$", "%d (%.1f%%)" %(diffs["bifidobacteriumAnimalis"][1][0], calcPerc(diffs["bifidobacteriumAnimalis"][1][0], iSecs["bifidobacteriumAnimalis"][0][0])), "Panaroo", \
		"%d (%.1f%%)" %(diffs["bifidobacteriumAnimalis"][1][1], calcPerc(diffs["bifidobacteriumAnimalis"][1][1], iSecs["bifidobacteriumAnimalis"][1][0])), "%d (%.1f%%)" %(diffs["yersiniaPestis"][1][0], \
		calcPerc(diffs["yersiniaPestis"][1][0], iSecs["yersiniaPestis"][0][0])), "Panaroo", "%d (%.1f%%)" %(diffs["yersiniaPestis"][1][1], calcPerc(diffs["yersiniaPestis"][1][1], \
		iSecs["yersiniaPestis"][1][0])), "%d (%.1f%%)" %(diffs["enterococcusFaecium"][1][0], calcPerc(diffs["enterococcusFaecium"][1][0], iSecs["enterococcusFaecium"][0][0])), "Panaroo", "%d (%.1f%%)" \
		%(diffs["enterococcusFaecium"][1][1], calcPerc(diffs["enterococcusFaecium"][1][1], iSecs["enterococcusFaecium"][1][0])), "%d (%.1f%%)" %(diffs["listeriaMonocytogenes"][1][0], calcPerc(\
		diffs["listeriaMonocytogenes"][1][0], iSecs["listeriaMonocytogenes"][0][0])), "Panaroo"]))
	print("\t".join(["", "%d (%.1f%%)" %(diffs["bifidobacteriumAnimalis"][2][0], calcPerc(diffs["bifidobacteriumAnimalis"][2][0], iSecs["bifidobacteriumAnimalis"][0][1])), "%d (%.1f%%)" \
		%(diffs["bifidobacteriumAnimalis"][2][1], calcPerc(diffs["bifidobacteriumAnimalis"][2][1], iSecs["bifidobacteriumAnimalis"][1][0])), "SibeliaZ", "%d (%.1f%%)" %(diffs["yersiniaPestis"][2][0], \
		calcPerc(diffs["yersiniaPestis"][2][0], iSecs["yersiniaPestis"][0][1])), "%d (%.1f%%)" %(diffs["yersiniaPestis"][2][1], calcPerc(diffs["yersiniaPestis"][2][1], iSecs["yersiniaPestis"][1][0])), \
		"SibeliaZ", "%d (%.1f%%)" %(diffs["enterococcusFaecium"][2][0], calcPerc(diffs["enterococcusFaecium"][2][0], iSecs["enterococcusFaecium"][0][1])), "%d (%.1f%%)" \
		%(diffs["enterococcusFaecium"][2][1], calcPerc(diffs["enterococcusFaecium"][2][1], iSecs["enterococcusFaecium"][1][0])), "SibeliaZ", "", ""]))
