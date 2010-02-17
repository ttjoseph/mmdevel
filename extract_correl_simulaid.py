#!/usr/bin/python
#
# Extracts the correlation matrix from Simulaid's correlation matrix output file.
# Starts with "The correlation matrix:" and ends with a blank line.

import sys
from math import sqrt

keep = False
data = []
for line in sys.stdin:
    if line.find("The correlation matrix") >= 0:
        keep = True
        continue
    elif len(line) <= 1:
	keep = False

    if keep: # format is 10i8
        for i in xrange(0, 80, 8):
            x = line[i:i+8].strip()
            if len(x) > 0:
	        data.append(float(x))

num_residues = int(sqrt(len(data)))
print >>sys.stderr, "%d numbers, which means %d residues." % (len(data), num_residues)
for i in xrange(num_residues):
    print " ".join([str(x) for x in data[i*num_residues:i*num_residues+num_residues]])
