#!/usr/bin/env python
#
# Generates a map of atoms from one PDB to another.
# The initial use case for this is modifying an existing PSF/PDB partway through a simulation
# so we can convert a regular lipid into a dual-topology one (or vice-versa, when lambda is
# 0 or 1). We need this sort of map so we can also modify the .coor and .vel files to continue
# the simulation without having to re-minimize and equilibrate and everything.
#
# Usage: <pdb1> <pdb2> <coorvel_prefix> <out_prefix>
# 
# pdb1 is the original PDB that you've been simulating, and for which you have good .coor and .vel
# files (specified by coorvel_prefix). pdb2 is the new PDB, that you've probably generated from 
# that .coor file after slight modification.
#
# So, the workflow you probably want is as follows:
#
# 1. Simulate some using original PSF/PDB, resulting in prod1.restart.coor and prod1.restart.vel.
# 2. Convert prod1.restart.coor into a PDB and edit it as you see fit. This probably involves
#    modifying the lipid head, and if this results in a charge state change, removing an ion
#    to maintain neutral charge in the system.
# 3. Use regenerate_psf.py to make a new PSF from that PDB.
# 4. Use this script, with the original and edited PDBs from prod1.restart.coor, to generate
#    updated .coor and .vel files.
# 5. Now you have everything you need to continue the simulation as modified.
#
# Tom Joseph, University of Pennsylvania
import sys
import argparse
import numpy as np
from libttj import CoorVel, ShadyPDB

def calculate_atoms_bijection(u1, u2):
    # print('calculate_atoms_bijection: hashing...', file=sys.stderr)
    u1_dict = hash_shadypdb(u1)
    u2_dict = hash_shadypdb(u2)
    # The number of keys in each dict should be the same as the number of atoms
    # If this is not true then we have duplicate atoms
    if len(u1_dict.keys()) != len(u1.atoms):
        print('Error: Degenerate atoms (with same name/coords) in first pdb', file=sys.stderr)
        sys.exit(1)
    if len(u2_dict.keys()) != len(u2.atoms):
        print('Error: Degenerate atoms (with same name/coords) in second pdb', file=sys.stderr)
        sys.exit(1)
    # print('calculate_atoms_bijection: done hashing', file=sys.stderr)

    a_to_b, a_not_in_b = calculate_atoms_injection(u1_dict, u2_dict)
    b_to_a, b_not_in_a = calculate_atoms_injection(u2_dict, u1_dict)

    # TODO: Sanity check: ensure the inverse of a_to_b == b_to_a
    # inverted_a_to_b = {v: k for k, v in a_to_b.items()}

    return a_to_b, np.array(a_not_in_b), np.array(b_not_in_a)


def hash_shadypdb(u):
    """Makes a dict out of a ShadyPDB, hashing atom name and position to the atom object itself."""
    d = {}
    for a_idx in range(len(u.atoms)):
        a = u.atoms[a_idx]
        d[f'{a.atomname} {a.x} {a.y} {a.z}'] = a_idx
    return d


def calculate_atoms_injection(u1_dict, u2_dict):
    a_to_b = {}
    a_not_in_b = []

    for a_hash, a_ix in u1_dict.items():
        if a_hash in u2_dict:
            a_to_b[a_ix] = u2_dict[a_hash]
        else:
            a_not_in_b.append(a_ix)
    # print('calculate_atoms_injection: done making map', file=sys.stderr)
    return a_to_b, a_not_in_b


def main():
    ap = argparse.ArgumentParser(description='Remap .coor and .vel files according to coorespondence between PDBs')
    ap.add_argument('pdb1', help='A PDB file')
    ap.add_argument('pdb2', help='The other PDB file. This is supposed to be a bijection')
    ap.add_argument('coorvel_prefix', help='.coor/.vel prefix')
    ap.add_argument('out_prefix', help='Out prefix for new .coor/.vel files that will work with pdb2')
    args = ap.parse_args()

    pdb1 = ShadyPDB(args.pdb1)
    pdb2 = ShadyPDB(args.pdb2)

    pdb1_to_pdb2, pdb1_only_atomidx, pdb2_only_atomidx = calculate_atoms_bijection(pdb1, pdb2)

    print(f'Atoms in pdb1 only: {pdb1_only_atomidx.size} of {len(pdb1.atoms)}')
    print(f'Atoms in pdb2 only: {pdb2_only_atomidx.size} of {len(pdb2.atoms)}')

    coor, vel = CoorVel(), CoorVel()
    coor.load(f'{args.coorvel_prefix}.coor')
    vel.load(f'{args.coorvel_prefix}.vel')

    out_coor, out_vel = CoorVel(), CoorVel()
    out_coor.blank(len(pdb2.atoms))
    out_vel.blank(len(pdb2.atoms))

    # Since pdb1_to_pdb2 is a bijection, we can just iterate to copy
    for a_ix, b_ix in pdb1_to_pdb2.items():
        a_offset, b_offset = a_ix*3, b_ix*3
        out_coor.coords[b_offset:b_offset+3] = coor.coords[a_offset:a_offset+3]
        out_vel.coords[b_offset:b_offset+3] = vel.coords[a_offset:a_offset+3]

    # TODO: fill in velocities we didn't copy
    # Interpolate, I guess? Or just give them a very low velocity and let the thermostat work its magic?

    # Tell user if they are adding atoms because we are leaving the velocities zero.
    if len(pdb2_only_atomidx) > 0:
        print(f'Warning: I am using {args.pdb2} for coordinates for its atom IDs:\n{pdb2_only_atomidx+1}', file=sys.stderr)
        print(f'Warning: And setting velocities to zero.', file=sys.stderr)
        print(f'Warning: This means you better have set the coordinates correctly in {args.pdb2}!', file=sys.stderr)
        if(len(pdb2_only_atomidx) > 20):
            print(f'Warning: Also, you are adding {len(pdb2_only_atomidx)} atoms, which is a lot of atoms to have zero velocities for.',
                file=sys.stderr)
        
    for b_ix in pdb2_only_atomidx:
        out_coor.coords[b_ix*3:b_ix*3+3] = pdb2.atoms[b_ix].position

    out_coor.write(f'{args.out_prefix}.coor')
    out_vel.write(f'{args.out_prefix}.vel')

    print(f"Wrote {out_coor.num_atoms} atoms to {args.out_prefix}.coor and {args.out_prefix}.vel.", file=sys.stderr)
    print(f"Ensure you generated a PSF to go with {args.pdb2}, perhaps using regenerate_psf.py,\nand then you ought to be able to run your simulation.",
        file=sys.stderr)


if __name__ == "__main__":
    main()
