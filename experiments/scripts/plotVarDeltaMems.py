#!/usr/bin/env python3

import argparse as args
from matplotlib import pyplot as plt

#This scripts plots the memory usages for varying delta values

def loadMemUsage(filename):
    for l in open(filename, 'r'):
        if l.find("Maximum") >= 0:
            return int(l.split(": ")[1])

if __name__ == '__main__':
    #Setting up the argument parser
    parser = args.ArgumentParser(description="This scripts plots the memory usages for varying delta values.")
    parser.add_argument('b', metavar='Benchmark', type=str, nargs='+', help="Benchmark files")
    parser.add_argument('-o', metavar='outFile', type=str, required=True, help="Output PDF name")

    arguments = parser.parse_args()
    mems = {}

    for f in arguments.b:
    	dlt = int(f.split('_d')[1].split('.txt')[0])
    	mems[dlt] = loadMemUsage(f)

    sortedDeltas = sorted(mems.keys())
    plt.style.use("grayscale")
    plt.plot(sortedDeltas, [mems[d] / 1000. for d in sortedDeltas], 'o-')
    plt.xlabel(r"$\delta$")
    plt.ylabel("Memory Peak (Mb)")
    plt.ylim([0, 850])
    plt.savefig(arguments.o, format="pdf")
