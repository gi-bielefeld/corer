#!/usr/bin/env python3

from Bio import SeqIO
import argparse as args
from Bio.SeqRecord import SeqRecord
from os.path import basename

#This script takes a genome assembly and gene annotations of this assembly and outputs a FASTA file with all gene sequences

#Setting up the argument parser
parser = args.ArgumentParser(description="This script creates a FASTA file containing all gene sequences for a given assembly")
parser.add_argument('-a', metavar='Assembly', type=str, required=True, help="Genome assembly in FASTA format")
parser.add_argument('-p', metavar='GenePreds', type=str, required=True, help="Gene predictions in GFF format")
parser.add_argument('-o', metavar='OutFile', type=str, required=True, help="Name of output file")
arguments = parser.parse_args()

#Load assembly
assembly = {c.id: c.seq for c in SeqIO.parse(arguments.a, "fasta")}
#Load gff file
isPredLine = False
records = []
counter = 0

for l in open(arguments.p):
    #Dicide how to parse the annotation file
    if l.find("AUGUSTUS") > -1:
        if not l.startswith('#'):
            elems = l.split('\t')

            if elems[2] == "gene":
                counter += 1
                records.append(SeqRecord(assembly[elems[0]][int(elems[3]) - 1:int(elems[4])], id=elems[0] + '_' + str(counter)))
    else:
        if l.startswith("##sequence-region"):
            isPredLine = True
            continue
        
        if l.startswith("##FASTA"):
            isPredLine = False
            continue
            
        if isPredLine:
            elems = l.split('\t')
            
            if elems[8].find("ID=") > -1:
                name = elems[8].split("ID=")[1].split(';')[0]
                start = int(elems[3]) - 1
                end = int(elems[4])
                # strand = elems[6]
                records.append(SeqRecord(assembly[list(assembly.keys())[0]][start:end], id=name))

# #Parse gene presents absence file
# frst = True
# coreGenes = []
# abs_qrm = ceil(arguments.q * 18)
# assemblyId = basename(arguments.a).split(".fasta")[0]

# for l in open(arguments.c, 'r'):
#     l = l.strip()
#     columns = l.split(',')
    
#     if frst:
#         colIdx = columns.index(assemblyId)
#         frst = False
#         continue
    
#     if columns[3:].count('') <= 18 - abs_qrm:
#         for g in columns[colIdx].split(';'):
#             if g != '' and g.find("refound") < 0:
#                 coreGenes.append(g)
                
# #Sanity check
# if len(coreGenes) != len(unique(coreGenes)):
#     print("WARNING: There is a gene occurring twice!", file=stderr)

# records = []

# for g in coreGenes:
#     subs = assembly.seq[predCoords[g][0]:predCoords[g][1]]
    
#     #Actually the strand doesn't matter...
#     #if predCoords[g][2] == '+':
#     #    seq = Seq(subs)
#     #else:
#     #    seq = Seq(Seq(subs).reverse_complement())

#     records.append(SeqRecord(subs, id=g))

#Write FASTA with gene sequences    
SeqIO.write(records, arguments.o, "fasta")