#!/usr/bin/env python
#
# How near does the specified ligand get to the specified other residues?
#
# Usage: ligand_near.py <ligand-spec> <residue-spec-1,residue-spec-2,...> <psf> <dcd> [dcd...]
#
# Where ligand-spec looks like HETA:1
# And residue-spec looks like PROA:182,PROB:182
import argparse
import sys
import numpy as np
from scipy.spatial.distance import euclidean
import MDAnalysis as mda

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Determines distance of ligand to residues). Outputs a CSV file for your plotting pleasure.')
    ap.add_argument('ligand_spec', help='Ligand residue ID, prefixed by "SEGID:".' )
    ap.add_argument('residues_spec', help='Binding site residue ID, prefixed by "SEGID:". To specify more than one, separate by commas, but no whitespace' )
    ap.add_argument('psf', help='PSF topology file')
    ap.add_argument('dcd', nargs='+', help='DCD trajectory file[s]')
    args = ap.parse_args()

    # Have to load the molecular system before we can select atoms
    print >>sys.stderr, 'Loading %s and %s...' % (args.psf, args.dcd)
    u = mda.Universe(args.psf, args.dcd)

    # Find the ligand
    ligand_segid, ligand_resid = args.ligand_spec.split(':')
    ligand_resid = int(ligand_resid)
    ligand = u.select_atoms('resid %d and segid %s' % (ligand_resid, ligand_segid))

    # Find each residue
    residue_strs = args.residues_spec.split(',')
    residues = []
    for spec in residue_strs:
        residue_segid, residue_resid = spec.split(':')
        residue_resid = int(residue_resid)
        residues.append(u.select_atoms('protein and name CA and segid %s and resid %d' % (residue_segid, residue_resid)))

    # Use the residues spec as the CSV header.
    print args.residues_spec
    # Quite convenient we made the user use commas eh?

    for ts in u.trajectory:
        vals = []
        for residue in residues:
            vals.append(euclidean(ligand.center_of_mass(), residue[0].position))
        print ','.join(['%.2f' % v for v in vals])