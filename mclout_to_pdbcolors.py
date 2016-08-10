#!/usr/bin/env python
# Converts the first line of mcl output into output suitable for
# input to colorize_pdb.py. Reads from stdout, writes to stdout.
import sys

# Get residue numbers that should be lit up.
# These are numbered starting at 1.
m = [int(i) for i in sys.stdin.readline().split()]

for i in xrange(max(m)):
    # Hilariously inefficient, but who cares?
    if i in m:
        print "1.0 ",
    else:
        print "0.0 ",
        