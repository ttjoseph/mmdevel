#!/usr/bin/env python
#
# Grimy script to weighted-sum the dA/dl from a NAMD simulation log, where
# colvars was used.
# Looks for lines containing both 'colvars:' and 'Lambda=', and extracts
# lambda and energy values.
import sys
from math import fabs

if __name__ == '__main__':
    rawlines = sys.stdin.readlines()
    vals = []
    # Save only the colvars TI lines
    for rawline in rawlines:
        if 'colvars:' in rawline and 'Lambda=' in rawline:
            (_, _, lambd, _, energy) = rawline.split()
            vals.append((float(lambd), float(energy)))
    # Assume energy is zero for the last point
    # vals.append((0.0, 0.0))

    total = 0.0
    for i in range(1, len(vals)):
        l0, e0 = vals[i-1]
        l1, e1 = vals[i]
        # Take the weight to be the difference between the two adjacent lambda values
        w = fabs(l1 - l0)
        # Take the average of the energies at the lambda points
        en = (e0 + e1) / 2.0
        print >>sys.stderr, 'Weight %.2f for energy %.2f kcal/mol at lambdas %.3f to %.3f' % (w, en, l0, l1)
        total += w * en

    print 'Total: %.2f kcal/mol' % total
