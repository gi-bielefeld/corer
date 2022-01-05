#!/usr/bin/env python3

import sys

#This script takes an annotation file of augustus and outputs a multiple fasta file containing amino acid sequences of all annotations

from os.path import basename
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

gname = ""
counter = 0
acc = basename(sys.argv[1]).split(".gff")[0]
isProtSeq = False
records = []

for l in open(sys.argv[1], 'r'):
    l = l.strip()
    
    if not l.startswith('#'):
        elems = l.split('\t')
       
        if elems[2] == "gene":
            gname = elems[0]
            counter += 1
            
    if l.startswith("# protein sequence = ["):
        isProtSeq = True
        seq = l.split('[')[1].split(']')[0]
    
        if l.find(']') > -1:
            records.append(SeqRecord(Seq(seq), id=acc + gname + '_' + str(counter)))
            isProtSeq = False
            
        continue
            
    if isProtSeq:
        seq += l.split(' ')[1].split(']')[0]
        
        if l.find(']') > -1:
            records.append(SeqRecord(Seq(seq), id=acc + gname + '_' + str(counter)))
            isProtSeq = False
            
SeqIO.write(records, sys.argv[1].split(".gff")[0] + ".fasta", "fasta")
