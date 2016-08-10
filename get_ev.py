#!/usr/bin/python
#
# Formats specified eigenvector with 3 coordinates per line (for use with ev_animate_pdb.py et al)
import sys
import os

if __name__ == '__main__':
    # Usage example (to get 2nd EV):
    # get_ev.py CA.PCA.eigenvectors.dat 2
    evnum = int(sys.argv[2])
    f = open(sys.argv[1])
    data = f.readlines()
    ev = [float(x) for x in data[evnum].strip().split()]
    for i in xrange(0, len(ev), 3):
        print "%f %f %f" % (ev[i], ev[i+1], ev[i+2])
