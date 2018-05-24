#!/usr/bin/env python
import argparse
import sys
import numpy as np
import MDAnalysis as mda

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Calculate distance between two residues over a trajectory')
    ap.add_argument('groups', help='Residue IDs, prefixed by "SEGID:", separated by hyphens, groups separated by commas, and no whitespace: e.g. "PROA:a-b-c" for the a-b-c angle.' )
    ap.add_argument('psf', help='PSF topology file')
    ap.add_argument('dcd', nargs='+', help='DCD trajectory file[s]')
    args = ap.parse_args()

    print >>sys.stderr, 'Loading %s and %s...' % (args.psf, args.dcd)
    u = mda.Universe(args.psf, args.dcd)
    print >>sys.stderr, 'Segments:', u.segments

    # Parse resid group list
    group_strs = args.groups.split(',')
    groups, segids = [], []
    for s in group_strs:
        segid, resid_str = s.split(':')
        segids.append(segid)
        resids = [int(a) for a in resid_str.split('-')]
        group = []
        for resid in resids:
            atoms = u.select_atoms('protein and segid %s and resid %d and name CA' % (segid, resid))
            if len(atoms) == 0:
                print >>sys.stderr, 'Error: Could not find the CA atom for segid %s and resid %d' % (segid, resid)
                exit(1)
            group.append(atoms[0])
        groups.append(group)

    # Convenient CSV header was provided by user
    print args.groups

    # Iterate through trajectory frames and calculate those angles
    for ts in u.trajectory:
        vals = []
        for group in groups:
            vals.append(np.linalg.norm(group[0].position - group[1].position))
        print ','.join(['%.2f' % v for v in vals])
