#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Computes how often any specified ligand was near the C-alpha atom of each residue.
# Spits out a histogram data file.
# This is useful for making sense of flooding simulations.

from __future__ import print_function
import argparse
import sys
from tqdm import tqdm
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import *


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Where do ligands go? Useful for flooding simulations')
    ap.add_argument('-c', '--cutoff', default=5, type=float, help='Distance cutoff in Angstroms')
    ap.add_argument('-t', '--percent-frames-threshold', default=0, type=float,
        help='Only print residues with this percent occupancy or greater (out of 100)')
    ap.add_argument('-n', '--num-frames', type=int,
        help='Limit the number of frames processed (beginning at the first frame)')
    ap.add_argument('-s', '--atomselect', action='store_true', default=False,
                    help='ligand_resname argument is actually an MDAnalysis atomselect string')
    ap.add_argument('-g', '--group-across-chains', action='store_true', default=False,
                    help='Group counts across chains (useful when the protein is multiple identical monomers)')
    ap.add_argument('--colorize-pdb-filename', help='Filename for output suitable for colorize_pdb.py (percent occupancy cutoff ignored)')
    ap.add_argument('ligand_resname', help='Residue name of ligand of interest (e.g. APM)')
    ap.add_argument('psf', help='PSF topology file')
    ap.add_argument('dcd', nargs='+', help='Trajectory file(s)')
    args = ap.parse_args()

    u = mda.Universe(args.psf, args.dcd)
    u.atoms.wrap(compound='fragments', inplace=True)
    protein_ca = u.select_atoms('protein and name CA')
    ligand = u.select_atoms(f'resname {args.ligand_resname}' if args.atomselect is False else args.ligand_resname)

    num_residues = len(protein_ca)
    first_resid = protein_ca.residues[0].resid
    num_frames = args.num_frames or len(u.trajectory)
    print(f'Atoms in ligand: {len(ligand)}', file=sys.stderr)
    print('Trajectory files: %s' % ', '.join(args.dcd), file=sys.stderr)
    print('Residues: %d' % num_residues, file=sys.stderr)
    print('First resid: %d' % first_resid, file=sys.stderr)
    print('Frames: To process %d of %d' % (num_frames, len(u.trajectory)), file=sys.stderr)
    print('Distance cutoff: %.2f Å for resname %s' % (args.cutoff, args.ligand_resname), file=sys.stderr)
    print('Percent occupancy threshold: %.2f%%' % args.percent_frames_threshold, file=sys.stderr)
    counts = np.zeros(num_residues)
    for ts in tqdm(u.trajectory, unit='frames'):
        if ts.frame >= num_frames:
            print(f"Stopping here since you only asked for the first {num_frames} frames.", file=sys.stderr)
            break

        d = distance_array(protein_ca.positions, ligand.positions, box=ts.dimensions)
        resid = 0
        for ca in d:
            min_dist = np.min(ca)
            if min_dist < args.cutoff:
                counts[resid] += 1
            resid += 1

    # We don't know beforehand how long the "coalesced" counts array needs to be,
    # and we need to index by resname+resid, so we use a dict.
    if args.group_across_chains is True:
        counts_coalesced = dict()
        for i in range(len(counts)):
            key = '%s%d' % (protein_ca[i].resname, protein_ca[i].resid)
            if key not in counts_coalesced:
                counts_coalesced[key] = 0
            counts_coalesced[key] += counts[i]

        for key in sorted(list(counts_coalesced.keys()), key=lambda res: int(res[3:])):
            percent_frames = 100*float(counts_coalesced[key])/num_frames
            if counts_coalesced[key] > 0 and percent_frames > args.percent_frames_threshold:
                print('%d\t%.2f%%\t%s' % (counts_coalesced[key], percent_frames, key))
    else:
        for i in range(len(counts)):
            percent_frames = 100*float(counts[i])/num_frames
            if counts[i] > 0 and percent_frames > args.percent_frames_threshold:
                print('%d\t%.2f%%\t%s%d' % (counts[i], percent_frames, protein_ca[i].resname, protein_ca[i].resid))

    if args.colorize_pdb_filename is not None:
        with open(args.colorize_pdb_filename, 'w') as f:
            print(' '.join([str(x) for x in counts/np.max(counts)]), file=f)

