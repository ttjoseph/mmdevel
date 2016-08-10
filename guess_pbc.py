#!/usr/bin/env python
import MDAnalysis as mda
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Guess NAMD periodic boundary conditions, for box cells only.')
parser.add_argument('psf', type=str, help='PSF file for the system to be simulated')
parser.add_argument('pdb', type=str, help='PDB file for the system to be simulated')
args = parser.parse_args()

u = mda.Universe(args.psf, args.pdb)
all_atoms = u.select_atoms('all')

smallest, largest = np.zeros(3), np.zeros(3)
total = np.zeros(3)

for a in all_atoms:
    smallest = np.minimum(smallest, a.position)
    largest = np.maximum(largest, a.position)

bounds = largest - smallest
center = (largest + smallest) / 2

a, b, c = np.ceil(bounds)
# For some reason we only care about the z-component of the origin.
# This is what CHARMM-GUI does.
zcen = center[2]

print("""# Periodic boundary conditions
cellBasisVector1 %.1f 0.0 0.0
cellBasisVector2 0.0 %.1f 0.0
cellBasisVector3 0.0 0.0 %.1f
cellOrigin 0.0 0.0 %.2f""" % (a, b, c, zcen))
