#!/usr/bin/env python3

import argparse
import pickle
from os.path import isfile
import MDAnalysis as mda
import MDAnalysis.lib.distances as mdadist
import numpy as np

# AutoDock atom type data lifted from Vina source code atom_constants.h
# generated from edited AD4_parameters.data using a script, 
# then covalent radius added from en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
# Columns:
#     name: radius, depth, solvation parameter, volume, covalent radius
atom_kind_data = { 
    "C": [2.00000,    0.15000,   -0.00143,   33.51030,   0.77], #  0
    "A": [2.00000,    0.15000,   -0.00052,   33.51030,   0.77], #  1
    "N": [1.75000,    0.16000,   -0.00162,   22.44930,   0.75], #  2
    "O": [1.60000,    0.20000,   -0.00251,   17.15730,   0.73], #  3
    "P": [2.10000,    0.20000,   -0.00110,   38.79240,   1.06], #  4
    "S": [2.00000,    0.20000,   -0.00214,   33.51030,   1.02], #  5
    "H": [1.00000,    0.02000,    0.00051,    0.00000,   0.37], #  6
    "F": [1.54500,    0.08000,   -0.00110,   15.44800,   0.71], #  7
    "I": [2.36000,    0.55000,   -0.00110,   55.05850,   1.33], #  8
    "NA": [1.75000,    0.16000,   -0.00162,   22.44930,   0.75], #  9
    "OA": [1.60000,    0.20000,   -0.00251,   17.15730,   0.73], # 10
    "SA": [2.00000,    0.20000,   -0.00214,   33.51030,   1.02], # 11
    "HD": [1.00000,    0.02000,    0.00051,    0.00000,   0.37], # 12
    "Mg": [0.65000,    0.87500,   -0.00110,    1.56000,   1.30], # 13
    "Mn": [0.65000,    0.87500,   -0.00110,    2.14000,   1.39], # 14
    "Zn": [0.74000,    0.55000,   -0.00110,    1.70000,   1.31], # 15
    "Ca": [0.99000,    0.55000,   -0.00110,    2.77000,   1.74], # 16
    "Fe": [0.65000,    0.01000,   -0.00110,    1.84000,   1.25], # 17
    "Cl": [2.04500,    0.27600,   -0.00110,   35.82350,   0.99], # 18
    "Br": [2.16500,    0.38900,   -0.00110,   42.56610,   1.14]  # 19
}

hydrophobic_atom_types = set(['C', 'F', 'Cl', 'Br', 'I'])
acceptor_atom_types = set(['N', 'O',])
donor_atom_types = set(['N', 'O',])

# TODO: This is stupid
RADIUS, DEPTH, SOLV_PARAM, VOLUME, COVALENT_RADIUS = 0, 1, 2, 3, 4

autodock_vdw_radii = {x: atom_kind_data[x][RADIUS] for x in atom_kind_data}

GAUSS1_WEIGHT = -0.0356
GAUSS2_WEIGHT = -0.00516
REPULSION_WEIGHT = 0.840
HYDROPHOBIC_WEIGHT = -0.0351
HBOND_WEIGHT = -0.587
NROT_WEIGHT = 0.585

def is_acceptor(a):
    return a.type in ['OA', 'NA']

def is_bonded_to_HD(a):
    for b in a.bonds:
        if b.partner(a).type == 'HD':
            return True
    return False

def atom_type_for_hbond(a):
    """Partial reimplementation of Vina's model::assign_types()"""
    donor_N_or_O = is_bonded_to_HD(a)
    acceptor = is_acceptor(a)

    if a.type[0] not in ['N', 'O']:
        return 'P'

    if acceptor and donor_N_or_O:
        return 'DA'
    elif acceptor:
        return 'A'
    elif donor_N_or_O:
        return 'D'
    else:
        return 'P'


def could_hbond(a1, a2):
    kind1 = atom_type_for_hbond(a1)
    kind2 = atom_type_for_hbond(a2)

    if (kind1, kind2) == ('A', 'D') or (kind1, kind2) == ('D', 'A') or (kind1, kind2) == ('DA', 'DA'):
        return True
    return False


def main():
    ap = argparse.ArgumentParser(description='Calculate Vina score for a static conformation')
    ap.add_argument('receptor', help='PDBQT file for receptor')
    ap.add_argument('ligand', help='PDBQT file for ligand, positioned where you want it relative to receptor')
    args = ap.parse_args()

    # Load receptor and ligand
    # This will take damn near forever due to guessing bonds
    ligand_u = mda.Universe(args.ligand, guess_bonds=True, vdwradii=autodock_vdw_radii)
    receptor_u = mda.Universe(args.receptor, guess_bonds=True, vdwradii=autodock_vdw_radii)

    # Vina score is described in: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3041641
    # We are only calculating the c_inter component

    # Calculate distance matrix
    dists = mdadist.distance_array(receptor_u.atoms.positions, ligand_u.atoms.positions)

    # Calculate d_ij, which is just the distance minus the vdW radii of the two atoms
    # I'm sure this is incredibly slow and there is some fast numpy array way to do this.
    # Probably by pre-assembling arrays of vdW radii and doing array subtractions.
    dij = dists.copy()
    hydrophobic_mask = np.zeros_like(dists)
    hbond_mask = np.zeros_like(dists)
    for r_i in range(len(receptor_u.atoms)):
        for l_i in range(len(ligand_u.atoms)):
            r_atom = receptor_u.atoms[r_i]
            l_atom = ligand_u.atoms[l_i]
            r_vdw = atom_kind_data[r_atom.type][RADIUS]
            l_vdw = atom_kind_data[l_atom.type][RADIUS]
            dij[r_i, l_i] -= r_vdw
            dij[r_i, l_i] -= l_vdw
            # If both atoms are hydrophobic, we care about their hydrophobic interaction score
            if r_atom.type in hydrophobic_atom_types and l_atom.type in hydrophobic_atom_types:
                hydrophobic_mask[r_i, l_i] = 1
            # Determine which atoms are involved in H-bonds
            if could_hbond(r_atom, l_atom):
                hbond_mask[r_i, l_i] = 1

    print('There are {} atom pairs in total'.format(len(receptor_u.atoms)*len(ligand_u.atoms)))
    print('{} atom pairs have hydrophobic interactions'.format(int(np.sum(hydrophobic_mask))))
    print('{} atom pairs can H-bond'.format(int(np.sum(hbond_mask))))

    # Calculate gauss1 component of score. Weight = -0.0356
    #   gauss1(d) = exp(-(d/0.5)^2)
    gauss1 = np.exp(-np.square(dij/0.5))

    # Calculate gauss2 component of score. Weight = -0.00516
    #   gauss2(d) = exp(-((d-3)/2)^2)
    gauss2 = np.exp(-np.square((dij-3.0)/2))

    # Calculate repulsion component of score. Weight = 0.840
    #   repulsion(d) = d^2 if d < 0, or 0 otherwise
    repulsion = np.square(np.clip(dij, a_min=None, a_max=0.0))

    # Calculate hydrophobic component of score. Weight = -0.0351
    # Both atoms must be hydrophobic or score is zero.
    #    hydrophobic(d) = 1 when d < 0.5
    #                   = 0 when d > 1.5
    #                   = linearly interpolated between these values otherwise
    hydrophobic = np.clip(1.5-dij, a_min=0.0, a_max=1.0)
    hydrophobic *= hydrophobic_mask

    # Calculate H-bonding component of score. Weight = -0.587
    #    hbond(d) = 1 when d < -0.7
    #             = 0 when d > 0
    #             = linearly interpolated between these values otherwise
    hbond = np.clip(-10*dij/7, a_min=0.0, a_max=1.0)
    hbond *= hbond_mask

    # Distance cutoff is 8A in Vina
    dist_cutoff_mask = np.zeros_like(dists)
    dist_cutoff_mask[dists <= 8] = 1

    gauss1 *= dist_cutoff_mask
    gauss2 *= dist_cutoff_mask
    repulsion *= dist_cutoff_mask
    hydrophobic *= dist_cutoff_mask
    hbond *= dist_cutoff_mask

    print('Gauss1: {}'.format(GAUSS1_WEIGHT*np.sum(gauss1)))
    print('Gauss2: {}'.format(GAUSS2_WEIGHT*np.sum(gauss2)))
    print('Repuls: {}'.format(REPULSION_WEIGHT*np.sum(repulsion)))
    print('Hphobc: {}'.format(HYDROPHOBIC_WEIGHT*np.sum(hydrophobic)))
    print('H-bond: {}'.format(HBOND_WEIGHT*np.sum(hbond)))

    c_inter = GAUSS1_WEIGHT*np.sum(gauss1) + \
        GAUSS2_WEIGHT*np.sum(gauss2) + \
        REPULSION_WEIGHT*np.sum(repulsion) + \
        HYDROPHOBIC_WEIGHT*np.sum(hydrophobic) + \
        HBOND_WEIGHT*np.sum(hbond)

    print(c_inter)

    # Calculate g(c_inter) as defined in the paper
    #    g(c_inter) = c_inter / (1 + 0.0585 * N_rot)
    #  where N_rot is the number of rotatable bonds between heavy atoms in the ligand

    # PLAN:
    #   Calculate difference between Vina scores in docked positions vs final MD positions for ligands
    #   embedded in lipid



if __name__ == '__main__':
    main()