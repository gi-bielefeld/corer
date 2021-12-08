#!/usr/bin/env python3

import numpy as np
import argparse as args
from matplotlib import pyplot as plt

#This script creates a peak memory plot from given benchmark files

SPECIES = ["bifidobacteriumAnimalis", "yersiniaPestis", "enterococcusFaecium", "listeriaMonocytogenes"]
TOOLS = ["Corer", "Panaroo", "SibeliaZ"]
LABELS = ["B. animalis", "Y. pestis", "E. faecium", "L. monocytogenes"]
BAR_WIDTH = 1

def loadMemUsage(filename):
    for l in open(filename, 'r'):
        if l.find("Maximum") >= 0:
            return int(l.split(": ")[1])

def mapBenchmarksToSpecies(speciesDict, benchmarkFiles):
    for f in benchmarkFiles:
        for s in SPECIES:
            if f.find(s) >= 0:
                speciesDict[s] = loadMemUsage(f)
                break

    return

if __name__ == '__main__':
    #Setting up the argument parser
    parser = args.ArgumentParser(description="This script creates a peak memory plot from given benchmark files.")
    parser.add_argument('-c', metavar='CorerMems', type=str, nargs='+', help="Benchmark files for Corer")
    parser.add_argument('-p', metavar='PanarooMems', type=str, nargs='+', help="Benchmark files for Panaroo")
    parser.add_argument('-s', metavar='SibeliaZMems', type=str, nargs='+', help="Benchmark files for SibeliaZ")
    parser.add_argument('-o', metavar='outFile', type=str, required=True, help="Output PDF name")

    arguments = parser.parse_args()

    toolbenchmarkFiles = [arguments.c, arguments.p, arguments.s]
    mem = {}

    for i in range(len(toolbenchmarkFiles)):
        if toolbenchmarkFiles[i]:
            mem[TOOLS[i]] = {}
            mapBenchmarksToSpecies(mem[TOOLS[i]], toolbenchmarkFiles[i])

    xCoords = np.array([len(mem) * (BAR_WIDTH + 0.25) * i for i in range(len(LABELS))])
    fig = plt.figure()
    plt.style.use("grayscale")

    for t in range(len(TOOLS)):
        if TOOLS[t] in mem:
            shift = t * BAR_WIDTH - ((len(TOOLS) / 2) - (BAR_WIDTH / 2))
            plt.bar(xCoords + shift, [mem[TOOLS[t]][p] for p in SPECIES if p in mem[TOOLS[t]]], BAR_WIDTH, label=TOOLS[t])
        
    plt.xticks(xCoords, LABELS)
    plt.ylabel("Memory Peak (Mb)")
    plt.legend()
    plt.savefig(arguments.o, format="pdf")