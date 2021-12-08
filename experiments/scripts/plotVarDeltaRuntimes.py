#!/usr/bin/env python3

import argparse as args
from matplotlib import pyplot as plt

#This scripts plots the given run times for varying delta values

def loadBenchmark(filename):
    for l in open(filename, 'r'):
        if l.find("User time") >= 0:
            return float(l.split(": ")[1])

if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This scripts plots the given run times for varying delta values.")
	parser.add_argument('b', metavar='Benchmark', type=str, nargs='+', help="Benchmark files")
	parser.add_argument('-o', metavar='outFile', type=str, required=True, help="Output PDF name")

	arguments = parser.parse_args()
	runtimes = {}

	for f in arguments.b:
		dlt = int(f.split('_d')[1].split('.txt')[0])
		runtimes[dlt] = loadBenchmark(f)

	sortedDeltas = sorted(runtimes.keys())
	plt.style.use("grayscale")
	plt.plot(sortedDeltas, [runtimes[d] for d in sortedDeltas], 'o-')
	plt.xlabel(r"$\delta$")
	plt.ylabel("User time (s)")
	plt.ylim([0,10000])
	plt.savefig(arguments.o, format="pdf")
