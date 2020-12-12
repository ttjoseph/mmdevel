#!/usr/bin/env python
import sys
import argparse
import csv
from pathlib import PurePath
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.dihedrals import Dihedral


def do_custom_selection(residue, names):
    sel = sum([residue.atoms[residue.atoms.names == n] for n in names])
    if len(sel) == 4:
        return sel
    else: # Invalid atom names...for this residue, but maybe not for others
        return None


# The MDAnalysis chi1_selection function is too stupid to understand all residues, so we roll our own
def do_chi1_selection(residue):
    sel = do_custom_selection(residue, ['N', 'CA', 'CB', 'CG'])
    if sel is None:
        # For Ile/Val
        sel = do_custom_selection(residue, ['N', 'CA', 'CB', 'CG1'])
    return sel


def main():
    ap = argparse.ArgumentParser(description='Calculate rotamer angles, defaulting to chi1 angles')
    ap.add_argument('selection', help='MDAnalysis atom selection expression for the residues you want, in quotes')
    ap.add_argument('psf', help='Topology for trajectory, probably in PSF format')
    ap.add_argument('dcd', nargs='+', help='Trajectory files, probably in DCD format')
    ap.add_argument('--custom-dihedral', help='A string with four atom names in the residue you want the dihedral of, such as "CA CB CG1 CD"')
    ap.add_argument('--skip-frames', default=0, type=int, help='Number of frames to skip (e.g. to take into account equilibration)')
    args = ap.parse_args()

    # Load and wrap trajectory, in case it was clipped across a sidechain
    u = mda.Universe(args.psf, args.dcd)
    u.atoms.wrap(compound='fragments', inplace=True)

    selection = u.select_atoms(args.selection)
    # We will cull this list of residues later, for when the custom residue spec only makes sense for some residues
    selected_residues = selection.residues
    # Chi1 dihedral angle is formed by N-CA-CB-CG
    if args.custom_dihedral is not None:
        custom_atomnames = args.custom_dihedral.strip().split()
        dihedral_sel = [do_custom_selection(residue, custom_atomnames) for residue in selected_residues]
    else:
        dihedral_sel = [do_chi1_selection(residue) for residue in selection.residues]

    # Remove null selections from this list of selections
    # Also culling the residue labels
    dihedral_sel_final = []
    selected_residues_final = []
    # print(dihedral_sel, file=sys.stderr)
    for i in range(len(dihedral_sel)):
        if dihedral_sel[i] is not None:
            dihedral_sel_final.append(dihedral_sel[i])
            selected_residues_final.append(selected_residues[i])

    # dihedral_sel = [x for x in dihedral_sel if x is not None]
    print(dihedral_sel_final, file=sys.stderr)
    dihedrals = Dihedral(dihedral_sel_final).run()

    labels = [f"{residue.resname}{residue.resid}" for residue in selected_residues_final]
    print(','.join(labels))
    writer = csv.writer(sys.stdout)
    for row in dihedrals.angles[args.skip_frames:]:
        writer.writerow(row)

if __name__ == '__main__':
    main()
