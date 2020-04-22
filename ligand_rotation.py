#!/usr/bin/env python
#
# Calculates the euler angles of the rotation of specified ligand over time.
# You could modify this to look at more than a ligand, I guess.
import csv
import sys
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from tqdm import tqdm
from scipy.spatial.transform import Rotation

def main():
    ap = argparse.ArgumentParser(description='Spit out CSV of ligand rotation Euler angles')
    ap.add_argument('ligand', help='Resname of ligand. Considers everything with this resname to be one single ligand, whether you like it or not.')
    ap.add_argument('psf', help='Topology file, probably in PSF format')
    ap.add_argument('dcd', nargs='+', help='Trajectory file(s), perhaps in DCD format')
    args = ap.parse_args()

    u = mda.Universe(args.psf, args.dcd)
    u.trajectory[0] # Select first frame
    ligand_u = u.select_atoms(f"resname {args.ligand}")
    ref = ligand_u.positions - ligand_u.center_of_mass()

    writer = csv.writer(sys.stdout)
    writer.writerow(['x', 'y', 'z'])

    for ts in tqdm(u.trajectory, unit='frames'):
        mobile0 = ligand_u.positions - ligand_u.center_of_mass()
        rot_mat, rmsd = align.rotation_matrix(mobile0, ref)
        rot = Rotation.from_matrix(rot_mat)
        writer.writerow(rot.as_euler('XYZ', degrees=True))



if __name__ == '__main__':
    main()