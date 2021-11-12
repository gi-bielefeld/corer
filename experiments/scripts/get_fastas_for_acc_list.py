#!/usr/bin/env python3

from sys import argv, stdout, exit
from Bio import SeqIO
from Bio import Entrez

Entrez.email='foo@bar.com'

def searchInDb(a_ID):
    handle = Entrez.efetch(db='nucleotide', id=a_ID, rettype='fasta')

    link = a_ID + ".fasta"
    local_file = open(link, 'w')
    local_file.write(handle.read())
    handle.close()
    local_file.close()

if __name__ == '__main__':
    if len(argv) != 2:
        print('\tusage: %s link to accession file missing')
        exit(1)
    name = argv[1]

    with open(name, "r") as ins:
    	for line in ins:
    		accession_ID = line.rstrip('\n')
    		print("Getting fasta file for", accession_ID)
    		searchInDb(accession_ID)