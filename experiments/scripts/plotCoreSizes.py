#!/usr/bin/env python3

import sys
from matplotlib import pyplot as plt

#This script plots the core sizes in k-mers with respect to the value of delta used

graphSizes = {}

for f in sys.argv[2:]:
	graphSizes = {int(f.split("_d")[1].split('.')[0]) = int(open(f, 'r').readline().split(' ')[10])}

keys = sorted(graphSizes.keys())
plt.style.use("grayscale")
plt.plot(keys, [graphSizes[k] / 1000. for k in keys], 'o-')
plt.xlabel(r"$\delta$")
plt.ylabel(r"$k$-mers")
plt.ylim([0, 12000])
plt.savefig(sys.argv[1], format="pdf")