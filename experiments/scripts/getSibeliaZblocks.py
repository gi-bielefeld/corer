#!/usr/bin/env python3

import argparse as args
from numpy import unique
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

#This script extracts blocks of sequences predicted by sibeliaZ from the corresponding assemblies

#Setting up the argument parser
parser = args.ArgumentParser(description="This script extracts blocks of sequences predicted by sibeliaZ from the corresponding assemblies.")
parser.add_argument('B', metavar='Blocks', type=str, help="Block coordinates file")
parser.add_argument('S', metavar='SeqFiles', type=str, nargs='+', help="Assemblies in FASTA format")
parser.add_argument('-q', metavar='quorum', type=int, required=True, help="Absolute quorum filter")
parser.add_argument('-o', metavar='outFile', type=str, default="blockSequences.fasta", help="Name of output file")

arguments = parser.parse_args()
genomeId = 0
genomeMapping = {}

#Read FASTA record names
for f in arguments.S:
	genomeId += 1

	for r in SeqIO.parse(f, "fasta"):
		genomeMapping[r.id] = [genomeId, r.seq]

#Read block coordinates
blocks = {}

for l in open(arguments.B, 'r'):
	if not l.startswith('#'):
		elems = l.split('\t')
		blockId = int(elems[8].split('=')[1])

		if not blockId in blocks:
			blocks[blockId] = []

		#Sanity check (Start position should never be larger than end position)
		if int(elems[4]) - int(elems[3]) < 0:
			print("ERROR: Start position larger than end position!", file=stderr)
			exit(-1)

		blocks[blockId].append([elems[0], int(elems[3]), int(elems[4]), elems[6]])

outrecords = []

for b in blocks:
	if arguments.q <= len(unique([genomeMapping[r[0]][0] for r in blocks[b]])):
		for r in blocks[b]:
			if r[3] == '+':
				outrecords.append(SeqRecord(genomeMapping[r[0]][1][r[1]:r[2]], id="%d:%s:%d:%d:%s" %(b, r[0], r[1], r[2], r[3])))
			else:
				#offsets seem to be wrong on the reverse complementary strand...
				start = r[1] - 1
				end = r[2] - 1
				outrecords.append(SeqRecord(genomeMapping[r[0]][1][start:end], id="%d:%s:%d:%d:%s" %(b, r[0], r[1], r[2], r[3])))

SeqIO.write(outrecords, arguments.o, "fasta")