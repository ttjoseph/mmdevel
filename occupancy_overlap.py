#!/usr/bin/env python
#
# Using the output of two ligand_occupancy.py runs, tells you which residues had the highest occupancy
# in both. This is useful for comparing two flooding simulations with the same protein but different
# ligands, such as with photoaffinity ligand and its native analog.
# Note: You probably want to use the '-g' option to ligand_occupancy.py so residues are coalesced
# across monomers.

import argparse
import sys
import operator
# from ligand_occupancy import residue_cmp

def occupancy_output_to_dict(lines):
    d = dict()
    for line in lines:
        count, percent, residue = line.split()
        if residue in d:
            print >>sys.stderr, 'Residue %s was mentioned twice, and I am not smart enough to deal with it' % residue
            exit(1)
        d[residue] = float(percent.strip('%'))
    return d

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Show overlapping residues from ligand_occupancy.py output')
    ap.add_argument('out1', help='First ligand_occupancy.py output file')
    ap.add_argument('out2', help='Second ligand_occupancy.py output file')
    args = ap.parse_args()

    with open(args.out1, 'r') as f:
        lines1 = f.readlines()

    with open(args.out2, 'r') as f:
        lines2 = f.readlines()

    # Convert each occupancy listing into a dict to make it easier to work with
    d1 = occupancy_output_to_dict(lines1)
    d2 = occupancy_output_to_dict(lines2)

    # Which residues show up in both places?
    intersecting_residues = list(set(d1.keys()) & set(d2.keys()))
    scores = dict()

    # We will be ranking residues by the product of their percentage occupancy
    for res in intersecting_residues:
        scores[res] = d1[res] * d2[res]

    # Sort by score descending, and print
    for res in sorted(scores.items(), key=operator.itemgetter(1), reverse=True):
        print '%s %.0f' % res