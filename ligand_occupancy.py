#!/usr/bin/env python
#
# Computes how often any specified ligand was near the C-alpha atom of each residue.
# Spits out a histogram data file.
# This is useful for making sense of flooding simulations.

import argparse
import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import *

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Where do ligands go? Useful for flooding simulations')
    ap.add_argument('-c', '--cutoff', default=5.0, type=float, help='Distance cutoff in Angstroms')
    ap.add_argument('ligand_resname', help='Residue name of ligand of interest (e.g. APM)')
    ap.add_argument('psf', help='PSF topology file')
    ap.add_argument('dcd', nargs='+', help='Trajectory file(s)')
    args = ap.parse_args()

    u = mda.Universe(args.psf, args.dcd)
    protein_ca = u.select_atoms('protein and name CA')
    ligand = u.select_atoms('resname %s' % args.ligand_resname)

    num_residues = len(protein_ca)
    first_resid = protein_ca.residues[0].resid
    print >>sys.stderr, 'Residues: %d' % num_residues
    print >>sys.stderr, 'First resid: %d' % first_resid
    counts = np.zeros(num_residues)
    for ts in u.trajectory:
        d = distance_array(protein_ca.positions, ligand.positions)
        resid = 0
        for ca in d:
            min_dist = np.min(ca)
            if min_dist < args.cutoff:
                counts[resid] += 1
            resid += 1

    for i in range(len(counts)):
        print '%d %d' % (first_resid+i, counts[i])
