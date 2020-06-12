#!/usr/bin/env python
#
# Automates the process of generating a flooded system which is hopefully suitable for feeding
# to CHARMM-GUI or similar. Useful for membrane proteins. You must start with a system
# previously generated with CHARMM-GUI, and plan to re-run its charmm scripts yourself.
#
# Requires packmol binary to be in your PATH.
import argparse
import sys
import subprocess
from tempfile import mkstemp
from os import unlink
import MDAnalysis as mda
import numpy as np


def main():
    ap = argparse.ArgumentParser(description='Generate a flood of ligands around a pre-oriented PDB')
    ap.add_argument('system_pdb', help='CHARMM-GUI step5_charmm2namd.pdb')
    ap.add_argument('ligand_pdb', help='PDB of one copy of ligand')
    ap.add_argument('--num_ligands', '-n', type=int, default=100, help='Number of ligand copies to start with')
    ap.add_argument('--around', type=float, default=7.0, help='How close a ligand can get to the solute')
    ap.add_argument('--out_pdb', default='flooded.pdb', help='Where to write the flooded ligands')
    ap.add_argument('--lipid_spec', default='resname POPC POPG POPE CHL1', help='Selection language for lipids (for MDAnalysis)')
    args = ap.parse_args()

    # Load the solute and see if there any DUM atoms, which OPM puts there to denote top and bottom of lipid bilayer
    system_u = mda.Universe(args.system_pdb)
    system_all = system_u.select_atoms('all')
    min_x, min_y, min_z = np.min(system_all.positions, axis=0)
    max_x, max_y, max_z = np.max(system_all.positions, axis=0)

    system_lipid = system_u.select_atoms(args.lipid_spec)
    if len(system_lipid) > 0:
        lower_membrane_z = np.min(system_lipid.positions[:,2])
        upper_membrane_z = np.max(system_lipid.positions[:,2])
        print(f'Detected lipid membrane with z bounds {lower_membrane_z:.2f} Å, {upper_membrane_z:.2f} Å',
              file=sys.stderr)
        exclude_bilayer_string = f'outside box {min_x} {min_y} {lower_membrane_z} {max_x} {max_y} {upper_membrane_z}'
    else:
        exclude_bilayer_string = ''

    # protein = system_u.select_atoms('protein and not resname DUM')
    ligand_u = mda.Universe(args.ligand_pdb)
    if len(ligand_u.residues) != 1:
        print(f'The ligand PDB {args.ligand_pdb} has {len(ligand_u.residues)} residues instead of just one', file=sys.stderr)
        print('I need it to be only one residue.')
        exit(1)

    ligand_resname = ligand_u.residues[0].resname
    ligands_already_there = system_u.select_atoms(f'resname {ligand_resname}')
    if len(ligands_already_there) > 0:
        print(f"Error: The system PDB {args.system_pdb} already has {ligand_resname} residues in it.", file=sys.stderr)
        print("Sorry, I don't know how to do a second flood.")
        exit(1)

    print(f'Flooding {args.system_pdb} with ligand {ligand_resname}.', file=sys.stderr)

    packmol_script = f"""tolerance 2.0
filetype pdb

output {args.out_pdb}

structure {args.ligand_pdb}
  number {args.num_ligands}
  inside box {min_x} {min_y} {min_z} {max_x} {max_y} {max_z}
  {exclude_bilayer_string}
end structure
"""
    _, tmpfile = mkstemp('.flood.packmol')
    with open(tmpfile, 'w') as f:
        f.write(packmol_script)

    # Run packmol to generate the ligand flood.
    # We have to actually use the < shell operator because packmol rewinds stdin,
    # and I guess the shell is smart enough to do that for redirected files but
    # not for true pipes.
    p = subprocess.run(f"packmol < {tmpfile}", shell=True, capture_output=True)
    unlink(tmpfile)
    if 'Success!' not in p.stdout.decode('utf-8'):
        print(p.stdout.decode('utf-8'), file=sys.stderr)
        print('generate-flood.py: Packmol failed in some way! See above.', file=sys.stderr)
        exit(1)

    # Combine the output of packmol with the solute
    flood_u = mda.Universe(args.out_pdb)
    both_u = mda.Merge(system_u.select_atoms('all'), flood_u.select_atoms('all'))

    # Exclude flooded ligands that clash with the protein
    good_ligands = both_u.select_atoms(f"same residue as (resname {ligand_resname} and not around {args.around} (protein or nucleic))")
    print('Left with', len(good_ligands.atoms),
          f'flooded ligand atoms ({len(good_ligands.residues)} residues) that are at least {args.around} Å from protein or nucleic.',
          file=sys.stderr)
    good_ligands.write(args.out_pdb)
    ligand_segid = ligand_u.residues[0].segid
    print(f"""
step1_pdbreader.inp
===================
You still have work to do. First, insert the following stanza in the correct spot in step1_pdbreader.inp:

! Flooded ligands
open read card unit 10 name {args.out_pdb}
read sequence pdb unit 10
generate {ligand_segid} setup warn first none last none
open read card unit 10 pdb name {args.out_pdb}
read coor pdb unit 10 resid

step3_size.inp (or equivalent)
==============================
There is a stanza in this script that determines the Z-axis extents of the protein, that would also now
count your new ligand flood and make the box too big. Edit so it looks something like:

! estimate protein extent along Z axis
coor stat sele .not. hydrogen .and. .not. resname {ligand_resname} end
calc protZmax = ?zmax
calc protZmin = ?zmin

toppar.str
==========
Finally, don't forget to also edit toppar.str so that it loads your RTF and PRM files for the ligand.
This may look something like:

open read card unit 10 name toppar/{ligand_resname.lower()}.rtf
read rtf card unit 10 append
open read card unit 10 name toppar/{ligand_resname.lower()}.prm
read para card unit 10 append flex

step5_charmm2namd.psf (or equivalent)
=====================================
This file isn't regenerated when you run all the CHARMM setup scripts, so you'll have to use
step5_assembly.psf instead it seems.
""")

if __name__ == '__main__':
    main()
