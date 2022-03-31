#!/usr/bin/env python3
#
# Renumbers protein residues in a PSF/PDB according to a template PDB.
#
# This is particularly useful when you've just used regenerate_psf.py, which uses
# VMD/psfgen for the heavy lifting, which doesn't preserve residue numbering.
#
# Tom Joseph, University of Pennsylvania
import argparse
import sys
from libttj import ShadyPSF, ShadyPDB, renumber_psf_pdb

def main():
    ap = argparse.ArgumentParser(description="Renumbers the solute in a NAMD PSF/PDB pair like a PDB template")
    ap.add_argument('psf', help='PSF file to renumber')
    ap.add_argument('pdb', help='PDB file to renumber')
    ap.add_argument('pdb_template', help='Template PDB')
    ap.add_argument('out_prefix', help='Prefix for output filenames')
    args = ap.parse_args()

    psf = ShadyPSF(args.psf)
    pdb = ShadyPDB(args.pdb)
    pdb_template = ShadyPDB(args.pdb_template)

    psf, pdb = renumber_psf_pdb(psf, pdb, pdb_template)

    # Write out newly modified PSF and PDB
    # For ShadyPDB, this is simply a matter of calling write_to_pdb()
    print(f'Info: Writing to {args.out_prefix}.psf and {args.out_prefix}.pdb', file=sys.stderr)
    psf.write_to_psf(f'{args.out_prefix}.psf')
    pdb.write_to_pdb(f'{args.out_prefix}.pdb')


if __name__ == "__main__":
    main()
 