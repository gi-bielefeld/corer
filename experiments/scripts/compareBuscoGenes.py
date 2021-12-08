#!/usr/bin/env python3

import argparse as args
import numpy as np

#This script compares the BUSCO genes found in the cores of assembly and read graphs

#Go through all genes that could be found in the core, check if they were mapped to a Busco id and collect all 
#Busco ids that could be found
def getCoreBuscos(queryRes, geneBuscoIdMap):
    coreBuscos = []

    for f in queryRes:
        acc = f.split('_q')[1].split(".tsv")[0]
        isFrst = True
    
        for l in open(f, 'r'):
            if isFrst:
                isFrst = False
                continue
            
            gene, flag = l.split('\t')
            geneKey = acc + gene
        
            if int(flag) > 0 and geneKey in geneBuscoIdMap:
                coreBuscos.append(geneBuscoIdMap[geneKey])
                
    return np.unique(coreBuscos)

#Setting up the argument parser
parser = args.ArgumentParser(description="This script compares the BUSCO genes found in the cores of assembly and read graphs.")
parser.add_argument('-b', metavar='BuscoGenes', type=str, required=True, nargs='+', help="Results of BUSCO gene analysis")
parser.add_argument('-a', metavar='AssemblyRes', type=str, required=True, nargs='+', help="Bifrost query results for assembly graph")
parser.add_argument('-r', metavar='ReadRes', type=str, required=True, nargs='+', help="Bifrost query results for read graph")

arguments = parser.parse_args()

#Go through all Busco result files and collect for each Busco id its corresponding gene predictions
geneBuscoIdMap = {}

for f in arguments.b:
    for l in open(f, 'r'):
        if not l.startswith('#'):
            elems = l.split('\t')
            
            #Check if mapping is "complete"
            if elems[1] == "Complete":
                geneBuscoIdMap[elems[2]] = elems[0]
       
assemblyCoreBuscos = getCoreBuscos(arguments.a, geneBuscoIdMap)
readCoreBuscos = getCoreBuscos(arguments.r, geneBuscoIdMap)

#How many BUSCO genes can be found among all predicted genes
print("The BUSCO gene analysis found %d BUSCO genes among all predictions" %len(np.unique([geneBuscoIdMap[g] for g in geneBuscoIdMap])))
#How many BUSCO genes were found in the assembly graph's core
print("%d BUSCO genes were found in the assembly graph's core" %len(assemblyCoreBuscos))
#How many BUSCO genes were found in the read graph's core
print("%d BUSCO genes were found in the read graph's core" %len(readCoreBuscos))
#How many BUSCO genes are found in the cores of both pangenomes?
print("%d BUSCO genes are found in both cores" %len(set(assemblyCoreBuscos).intersection(set(readCoreBuscos))))
