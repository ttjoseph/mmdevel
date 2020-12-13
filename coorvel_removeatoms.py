#!/usr/bin/env python
#
# Remove selected coordinates from NAMD binary .coor and .vel files.
# Makes the very significant assumption that everything is little endian.
#
# Tom Joseph, University of Pennsylvania
import sys
import argparse
import numpy as np
# from scipy.io import FortranFile


class CoorVel(object):
    """NAMD-format binary .coor or .vel file.

    These have the same format: number of atoms (int32), followed by xyz coordinates for each atom as doubles.
    Importantly, we assume everything is little endian, because most machines running NAMD these days are
    in fact little endian.
    """
    def __init__(self, fname):
        self.load(fname)

    def load(self, fname):
        """Loads .coor or .vel file into this object."""
        # So format is totatoms (32-bit integer) then 3 doubles for each atom
        # Assume little endian. Not sure if NAMD etc always write little endian or what.
        num_atoms = np.fromfile(fname, dtype=np.dtype('<i4'), count=1, offset=0)
        self.num_atoms = num_atoms[0]
        self.coords = np.fromfile(fname, dtype=np.dtype('<f8'), offset=4, count=self.num_atoms*3)

    def write(self, fname):
        """Writes this object in binary format suitable for NAMD."""
        f = open(fname, 'wb')
        f.write(self.num_atoms.astype(np.dtype('<i4')).tobytes())
        f.write(self.coords.tobytes())
        f.close()

    def excise(self, atom_indices):
        """Removes specified atoms in place."""
        # Convert atom indices to coordinate indices
        # Ex: [2, 3] -> [2*3, 2*3+1, 2*3+2, 3*3, 3*3+1, 3*3+2]
        x = atom_indices * 3
        y = atom_indices * 3 + 1
        z = atom_indices * 3 + 2
        xyz = np.concatenate((x, y, z), axis=None)
        xyz.sort()
        # print(f'Deleting coord indices: {xyz}', file=sys.stderr)
        # print(self.coords[xyz].reshape((-1, 3)), file=sys.stderr)
        self.coords = np.delete(self.coords, xyz)
        self.num_atoms -= len(atom_indices)



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

    coor = CoorVel(f'{args.prefix}.coor')
    vel = CoorVel(f'{args.prefix}.vel')
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