#!/usr/bin/env python
# Reads a square matrix of numbers on stdin and outputs an mcl --abc format
# file on stdout.
import sys

# The data itself!
m = []

# Read first line to see how big the matrix is supposed to be
tokens = sys.stdin.readline().split()
dimension = len(tokens)
m.extend(tokens)

# Read len(tokens)-1 lines
for i in xrange(dimension-1):
    m.extend(sys.stdin.readline().split())

# Convert to floating-point    
m = [float(i) for i in m]

# Addressing:
# m[dimension * row + col]

# Find all nonzero cells and dump them in mcl's format, e.g.:
# 23 518 0.2481
# We don't assume the matrix is symmetric, even though for the
# correlation analysis stuff it probably is
for row in xrange(dimension):
    for col in xrange(dimension):
        val = m[dimension * row + col]
        if val != 0: print "%d %d %f" % (row+1, col+1, val)