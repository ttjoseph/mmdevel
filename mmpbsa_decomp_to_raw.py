#!/usr/bin/env python
# Converts MMPBSA.py idecomp output to a raw text matrix suitable for loading into Matlab etc.
import sys

read_stdev = False

if len(sys.argv) > 2 and sys.argv[2] == 'stdev':
     print >>sys.stderr, "# Outputting standard deviation matrix!!"
     read_stdev = True

# Load file and strip whitespace and empty lines
lines = open(sys.argv[1]).readlines()
lines = [l.strip() for l in lines]
lines = [l for l in lines if len(l) > 0]
tmp = lines[-1].split()
numres1 = int(tmp[1])
numres2 = int(tmp[4])
print >>sys.stderr, "This file %s has a %d by %d matrix." % (sys.argv[1], numres1, numres2)
if numres1 != numres2:
    print >>sys.stderr, "I demand a square matrix!"
    sys.exit(1)
    
M = [float('nan')] * numres1 * numres2
for l in lines:
    try:
        tmp = l.split()
        i = int(tmp[1])-1
        j = int(tmp[4])-1
        total_energy = float(tmp[26])
        if read_stdev: total_energy = float(tmp[28])

        M[numres1*i+j] = total_energy
    except ValueError:
        # Error in conversion means this wasn't a valid line, so try again with the next line.
        # This is quick and dirty indeed.
        continue
    except TypeError:
        continue
    except IndexError:
        continue
        
for i in xrange(numres1):
    for j in xrange(numres2):
        if j != 0: print " ",
        print M[i*numres1+j],
    print ""
