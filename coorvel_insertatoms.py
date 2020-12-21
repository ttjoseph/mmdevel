#!/usr/bin/env python
#
# Inserts coordinates into .coor and .vel files.
# You probably want to do this if you are say adding an ion to a pre-existing simulation.
# Makes the very significant assumption that everything is little endian.
#
# Tom Joseph, University of Pennsylvania
import sys
import argparse
import numpy as np
import MDAnalysis as mda
from libttj import CoorVel


def main():
    ap = argparse.ArgumentParser(description='Insert atom coordinates from a PDB and guessed velocities into .coor/.vel')
    ap.add_argument('prefix', help='Prefix to .coor and .vel filenames')
    ap.add_argument('pdb', help='PDB containing coordinates to insert. Of course all else in the PDB is ignored')
    ap.add_argument('insert_after_atomid', type=int, help='Insert those atoms after this 1-based starting atom ID')
    ap.add_argument('out_prefix', help='Prefix for output .coor and .vel filenames')
    args = ap.parse_args()

    print(f'System endianness is {sys.byteorder}.', file=sys.stderr)
    if sys.byteorder != 'little':
        print(f'Error: But this script expects little-endian data. If you know what you\'re doing, edit this script and try again.')
        exit(1)

    coor = CoorVel(f'{args.prefix}.coor')
    vel = CoorVel(f'{args.prefix}.vel')
    print(f'Starting with {coor.num_atoms} atoms in {args.prefix}.coor and {vel.num_atoms} atoms in {args.prefix}.vel.', file=sys.stderr)

    if coor.num_atoms != vel.num_atoms:
        print(f'Error: {args.prefix}.coor contains {coor.num_atoms} atoms but {args.prefix}.vel contains {vel.num_atoms} atoms.', file=sys.stderr)
        print(f'Error: That is no good and suggests there is something sinister going on.', file=sys.stderr)
        exit(1)

    u = mda.Universe(args.pdb)
    print(f'You want to insert {len(u.atoms)} atom(s) after atom ID {args.insert_after_atomid} (1-based counting).')
    print(f'This will change the numbering, but atom ID {args.insert_after_atomid} will still point to the same atom.')
    # User gave us a 1-based index but CoorVel.insert expects a 0-based index.
    coor.insert(u.atoms, args.insert_after_atomid-1)
    # TODO: We can't just insert crap into the velocities. Probably just copy the one before it? Or just zeros?
    # vel.insert(u.atoms, args.insert_after_atomid-1)

    coor.write(f'{args.out_prefix}.coor')
    vel.write(f'{args.out_prefix}.vel')

    print(f"Wrote {coor.num_atoms} atoms to {args.out_prefix}.coor and {args.out_prefix}.vel.", file=sys.stderr)



if __name__ == "__main__":
    main()