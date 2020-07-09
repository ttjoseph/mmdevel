#!/usr/bin/env python
import sys
import argparse
import sys
import csv
from pathlib import PurePath
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.dihedrals import Dihedral


def main():
    ap = argparse.ArgumentParser(description='Analyze rotamers')
    ap.add_argument('selection', help='MDAnalysis atom selection expression for the residues you want, in quotes')
    ap.add_argument('psf', help='Topology for trajectory, probably in PSF format')
    ap.add_argument('dcd', nargs='+', help='Trajectory files, probably in DCD format')
    ap.add_argument('--skip-frames', default=0, type=int, help='Number of frames to skip (e.g. to take into account equilibration)')
    args = ap.parse_args()

    # Load and align trajectory to itself
    u = mda.Universe(args.psf, args.dcd)
    u.atoms.wrap(compound='fragments', inplace=True)
    # alignment = align.AlignTraj(u, u, in_memory=True)

    selection = u.select_atoms(args.selection)
    chi1s = [residue.chi1_selection() for residue in selection.residues]
    dihedrals = Dihedral(chi1s).run()

    labels = [f"{residue.resname}{residue.resid}" for residue in selection.residues]
    print(','.join(labels))
    writer = csv.writer(sys.stdout)
    for row in dihedrals.angles[args.skip_frames:]:
        writer.writerow(row)

if __name__ == '__main__':
    main()
