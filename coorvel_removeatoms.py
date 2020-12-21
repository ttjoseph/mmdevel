#!/usr/bin/env python
#
# Remove selected coordinates from NAMD binary .coor and .vel files.
# Makes the very significant assumption that everything is little endian.
#
# Tom Joseph, University of Pennsylvania
import sys
import argparse
import numpy as np
from libttj import CoorVel

def main():
    ap = argparse.ArgumentParser(description='Excise selected atoms from NAMD binary .coor/.vel files')
    ap.add_argument('prefix', help='Prefix to .coor and .vel filenames')
    ap.add_argument('start_atomid', type=int, help='1-based starting atom ID to include in excision')
    ap.add_argument('end_atomid', type=int, help='1-based ending atom ID to include in excision')
    ap.add_argument('out_prefix', help='Output prefix to name output .coor and .vel files')
    args = ap.parse_args()

    print(f'System endianness is {sys.byteorder}.', file=sys.stderr)
    if sys.byteorder != 'little':
        print(f'Error: But this script expects little-endian data. If you know what you\'re doing, edit this script and try again.')
        exit(1)

    coor = CoorVel()
    coor.load(f'{args.prefix}.coor')
    vel = CoorVel()
    vel.load(f'{args.prefix}.vel')
    print(f'Starting with {coor.num_atoms} atoms in {args.prefix}.coor and {vel.num_atoms} atoms in {args.prefix}.vel.', file=sys.stderr)

    if coor.num_atoms != vel.num_atoms:
        print(f'Error: {args.prefix}.coor contains {coor.num_atoms} atoms but {args.prefix}.vel contains {vel.num_atoms} atoms.', file=sys.stderr)
        print(f'Error: That is no good and suggests there is something sinister going on.', file=sys.stderr)
        exit(1)

    # print(f'Shape of coor.coords: {coor.coords.shape}', file=sys.stderr)
    # print(f'Shape of vel.coords: {vel.coords.shape}', file=sys.stderr)

    # Note this is an inclusive range that is 1-based. So we convert to 0-based
    # (and note arange is a half-interval)
    atom_indices = np.arange(args.start_atomid-1, args.end_atomid)
    print(f"Excising {len(atom_indices)} atoms with these 1-based IDs:\n{', '.join([str(x) for x in atom_indices.tolist()])}", file=sys.stderr)
    coor.excise(atom_indices)
    vel.excise(atom_indices)

    # print(f'Shape of coor.coords: {coor.coords.shape}', file=sys.stderr)
    # print(f'Shape of vel.coords: {vel.coords.shape}', file=sys.stderr)

    coor.write(f'{args.out_prefix}.coor')
    vel.write(f'{args.out_prefix}.vel')

    print(f"Wrote {coor.num_atoms} atoms to {args.out_prefix}.coor and {args.out_prefix}.vel.", file=sys.stderr)


if __name__ == "__main__":
    main()