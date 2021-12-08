#!/usr/bin/env python3

import sys

#This script takes blast result files and outputs for each query its maximum relative coverage by any alignment found

for f in sys.argv[1:]:
	currentQuery = ""
	lastWasQueryName = False
	maxAlgnLen = 0.0

	for l in open(f, 'r'):
	    l = l.strip()
	    
	    if l.startswith("Query= "):
	        if currentQuery != "":
	            print(currentQuery, maxAlgnLen / qLen)

	            maxAlgnLen = 0.0
	            
	        currentQuery = l.split(' ')[1]
	        lastWasQueryName = True
	        
	    if l.startswith("Length=") and lastWasQueryName:
	        lastWasQueryName = False
	        qLen = float(l.split('=')[1])
	        
	    if l.find("Identities") > -1:
	        maxAlgnLen = max(maxAlgnLen, float(l.split('/')[1].split(' ')[0]))

	print(currentQuery, maxAlgnLen / qLen)