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
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import re

# Sort key for a natural sort, so that "prod10" comes after "prod2"
# Code from: https://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower() for text in _nsre.split(s)]


def abbrev_to_full_atomspec(spec):
    tokens = spec.split(':')
    if len(tokens) != 2:
        return None
    segid = tokens[0]
    # Get atom names if any
    tokens = tokens[1].split('/')
    resid = tokens[0]
    if len(tokens) == 1:
        return ['segid {} and resid {}'.format(segid, resid),]
    
    atomnames = tokens[1].split('.')
    out_specs = []
    for a in atomnames:
        out_specs.append('segid {} and resid {} and name {}'.format(segid, resid, a))
    return out_specs


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Determines the nearest-atom distance of ligand to residues. Outputs a CSV file for your plotting pleasure.')
    ap.add_argument('ligand_spec', help='Ligand residue ID, prefixed by "SEGID:".' )
    ap.add_argument('residues_spec', help='Binding site residue ID, prefixed by "SEGID:". To specify more than one, separate by commas, but no whitespace' )
    ap.add_argument('psf', help='PSF topology file')
    ap.add_argument('dcd', nargs='+', help='DCD trajectory file[s]')
    args = ap.parse_args()

    # Have to load the molecular system before we can select atoms
    print >>sys.stderr, 'Loading %s and %s...' % (args.psf, ', '.join(sorted(args.dcd, key=natural_sort_key)))
    u = mda.Universe(args.psf, sorted(args.dcd, key=natural_sort_key))

    # Find the ligand
    #ligand_segid, ligand_resid = args.ligand_spec.split(':')
    #ligand_resid = int(ligand_resid)
    out_specs = abbrev_to_full_atomspec(args.ligand_spec)
    #ligand = u.select_atoms('resid %d and segid %s' % (ligand_resid, ligand_segid))
    if len(out_specs) != 1:
        print >>sys.stderr, 'Ligand spec is too ambiguous'
        exit(1)
    ligand = u.select_atoms(out_specs[0])

    # Find each residue
    residue_strs = args.residues_spec.split(',')
    field_names = []
    residues = []
    for spec in residue_strs:
        #residue_segid, residue_resid = spec.split(':')
        #residue_resid = int(residue_resid)
        #residues.append(u.select_atoms('protein and segid %s and resid %d' % (residue_segid, residue_resid)))
        out_specs = abbrev_to_full_atomspec(spec)
        for s in out_specs:
            print >>sys.stderr, '{} = {}'.format(spec, s)
            field_names.append(s)
            residues.append(u.select_atoms(s))

    # Use the residues spec as the CSV header.
    print ','.join(field_names)
    # Quite convenient we made the user use commas eh?

    for ts in u.trajectory:
        vals = []
        for residue in residues:
            da = distance_array(ligand.positions, residue.positions)
            vals.append(np.min(np.min(da)))
            # vals.append(euclidean(ligand.center_of_mass(), residue[0].position))
        print ','.join(['%.2f' % v for v in vals])