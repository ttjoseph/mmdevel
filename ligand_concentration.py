#!/usr/bin/env python
#
# Returns the concentration of a specified ligand over time.
# There can be only one ligand in the simulation.
import MDAnalysis as mda
import argparse

AVOGADRO = 6.022140857e23
CUBIC_ANGSTROMS_PER_LITER = 1e27

def main():
    ap = argparse.ArgumentParser(description='Return the concentration of a specified ligand over time.')
    ap.add_argument('ligand', help='Selection spec of the ligand in question')
    ap.add_argument('topology', help='Topology file of the system, probably in PSF format')
    ap.add_argument('trajectory', nargs='+', help='Trajectory file(s) to be analyzed')
    args = ap.parse_args()

    u = mda.Universe(args.topology, args.trajectory)    
    ligand = u.select_atoms(args.ligand)
    print('# Ligand is:')
    print(ligand)
    print('# Dimensions of box:')
    print(u.dimensions)

    if u.dimensions[3] != 90 or u.dimensions[4] != 90 or u.dimensions[5] != 90:
        print('Error: Simulation box is not rectangular. I am not smart enough to deal with that.')
        return 1

    for ts in u.trajectory:
        volume = u.dimensions[0] * u.dimensions[1] * u.dimensions[2]
        volumes_per_liter = CUBIC_ANGSTROMS_PER_LITER / volume
        # Since we assume one ligand per simulation volume, there must be volumes_per_liter
        # ligands in a liter's worth of this system. That, divided by Avogadro's number, is
        # the number of moles of the ligand in our notional liter.
        molar_concentration = volumes_per_liter / AVOGADRO
        print('%d: Volume: %.1f A^3 Conc: %.5f M' % (ts.frame, volume, molar_concentration))
        

if __name__ == '__main__':
    import sys
    sys.exit(int(main() or 0))