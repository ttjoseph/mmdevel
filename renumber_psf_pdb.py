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
from libttj import ShadyPSF, ShadyPDB

def main():
    ap = argparse.ArgumentParser(description="Renumbers the protein in a NAMD PSF/PDB pair like a PDB template")
    ap.add_argument('psf', help='PSF file to renumber')
    ap.add_argument('pdb', help='PDB file to renumber')
    ap.add_argument('pdb_template', help='Template PDB')
    args = ap.parse_args()

    psf = ShadyPSF(args.psf)
    pdb = ShadyPDB(args.pdb)
    pdb_template = ShadyPDB(args.pdb_template)

    # First, determine correspondence between pdb and pdb_template by coordinates, out to perhaps two decimal
    # places? If two atoms differ only by 0.009 Angstroms there is likely something wrong
    def coord_key(a):
        return f'{a.x:.2f} {a.y:.2f} {a.z:.2f}'
    
    # Invert atom -> coord maps
    coord_to_atom_pdb = {coord_key(a): a for a in pdb.atoms}
    coord_to_atom_pdb_template = {coord_key(a): a for a in pdb_template.atoms}

    if len(coord_to_atom_pdb) != len(pdb.atoms):
        print(f"Error: Only {len(coord_to_atom_pdb)} distinct coordinates in {args.pdb}, which means atoms overlap!",
            file=sys.stderr)
        sys.exit(1)
    if len(coord_to_atom_pdb_template) != len(pdb_template.atoms):
        print(f"Error: Only {len(coord_to_atom_pdb_template)} distinct coordinates in {args.pdb_template}, which means atoms overlap!",
            file=sys.stderr)
        sys.exit(1)

    num_atoms_not_found = 0
    for k, template_atom in coord_to_atom_pdb_template.items():
        if k in coord_to_atom_pdb:
            target_atom = coord_to_atom_pdb[k]
            psf_target_atom = psf.atoms[target_atom.atomindex]
            if psf_target_atom[ShadyPSF.RESNAME] != target_atom.resname or \
                psf_target_atom[ShadyPSF.RESID] != target_atom.resid:
                print(f'Error: Target PSF and PDB do not appear to correspond to the same structure!',
                    file=sys.stderr)
                sys.exit(1)

            # Don't renumber water
            if target_atom.resname in ['TIP3', 'TIP4', 'WAT', 'HOH']:
                continue

            # print(f'{target_atom.resname} {psf_target_atom[ShadyPSF.RESNAME]} {target_atom.resid} -> {template_atom.resname} {template_atom.resid}')
            # Now we need to locate this same atom in the target PSF.
            # We must assume the PDB and PSF have the same atoms in the same order.
            psf.atoms[target_atom.atomindex][ShadyPSF.RESID] = template_atom.resid
            pdb.atoms[target_atom.atomindex].resid = template_atom.resid

        else:
            num_atoms_not_found += 1

    print(f'Info: {num_atoms_not_found} atoms in template not found in target PDB', file=sys.stderr)
    # TODO: Write out newly modified PSF and PDB
    # For ShadyPDB, this is simply a matter of calling write_to_pdb()
    print(f'DEBUG: Writing to foo.psf and foo.pdb', file=sys.stderr)
    psf.write_to_psf('foo.psf')
    pdb.write_to_pdb('foo.pdb')


if __name__ == "__main__":
    main()
 