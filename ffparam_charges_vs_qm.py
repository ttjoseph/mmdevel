#!/usr/bin/env python
import sys
import argparse
import re
import yaml
from psf_prm_to_bin import get_psf_block_header, read_int_block
import csv

# Constants lifted from FFTK
HARTREES_PER_KCAL = 1.041308e-21
THINGS_PER_MOLE = 6.02214e23
COULOMB = 332.0636
# QM_ENERGY_SCALE = 1.16 # From FFTK, and currently a mystery
# QM_ENERGY_SCALE = (-6.95/-4.84) # From isoflurane paper


def parse_waterint_log(filename):
    """Extracts SCF energies from a water interaction calculation done by Gaussian."""
    lines = open(filename).readlines()
    energies, coords, optimized_params = [], [], []
    current_frame_coords = []
    current_optimized_params = {}
    optimized = False
    i = 0
    while i < len(lines):
        if re.search(r'Error termination', lines[i]):
            sys.stderr.write('Careful, %s has a Gaussian error termination.' % filename)
            return energies
        
        if re.search(r'Normal termination of Gaussian', lines[i]):
            # All is well
            # sys.stderr.write('Gaussian did not complain in %s. Nice!\n' % filename)
            break
            
        # Extract parameters that were optimized, such as distances and angles
        if re.search(r'Optimized Parameters', lines[i]):
            i += 5
            current_optimized_params = {}
            optimized = True
            while lines[i][1] != '-':
                tokens = lines[i].split()
                current_optimized_params[tokens[1]] = float(tokens[2])
                i += 1
        
        if re.search(r'Standard orientation:', lines[i]) or re.search(r'Z-Matrix orientation:', lines[i]):
            # Save the input coordinates for this piece of the scan
            # Eat the header
            i += 5
            # Read the coordinates
            current_frame_coords = []
            while lines[i][1] != '-':
                (idx, atomic_num, atomic_type, x, y, z) = lines[i].split()
                x, y, z = float(x), float(y), float(z)
                # Ignore dummy atoms which Gaussian says have atomic number -1
                if int(atomic_num) > 0:
                    current_frame_coords.append((int(atomic_num), x, y, z))
                i += 1

        # Look for the calculated QM energies
        scfdone_match = re.search(r'SCF Done:(.+)', lines[i])
        energy_match = re.search(r'Energy=(.+) NIter=(.+)', lines[i])
        scf = 0
        if scfdone_match or energy_match:
            tokens = lines[i].split()
            if scfdone_match:
                scf = float(tokens[4]) # e.g. E(RB3LYP) = <energy>
            else:
                scf = float(tokens[1])

        scf = scf * HARTREES_PER_KCAL * THINGS_PER_MOLE
        if scf != 0:
            # Keep the latest coordinate set for which there is an energy
            coords.append(current_frame_coords)
            if optimized is False:
                current_optimized_params = {}
            optimized_params.append(current_optimized_params)
            energies.append(scf)
            optimized = False
            # print '%s: Energy: %f (relative: %f)' % (filename, scf, scf - energies[0])

        i += 1
    print >>sys.stderr, 'Found %d energies and %d coord sets in %s' % (len(energies), len(coords), filename)
    return energies, coords, optimized_params


def load_psf(filename):
    psf = open(filename, "r")

    # Eat header
    l = psf.readline()
    if l[0:3] != "PSF":
        print >>sys.stderr, "%s doesn't look like a PSF file to me." % filename
        sys.exit(1)
    psf.readline() # Blank line
    (numlines, kind) = get_psf_block_header(psf)
    for i in xrange(numlines): psf.readline()
    psf.readline()

    # !NATOM block contains, for each atom, its type and charge, as well as residue assignment
    (num_atoms, kind) = get_psf_block_header(psf)
    print >>sys.stderr, "Number of atoms: %d" % num_atoms
    residue_map, charges, atom_types = [], [], []

    # Parse each atom line
    #          1         2         3         4         5         6
    #0123456789012345678901234567890123456789012345678901234567890123456789
    #       1 LIGH 1    ALA  N    NH3   -0.300000       14.0070           0
    found_end_of_solute = False
    num_solute_atoms = 0
    last_resid, last_segment = -1, None
    num_residues = num_solute_residues = 0
    cys_sg_list = [] # atom indices of CYS SG atoms
    for index in xrange(num_atoms):
        l = psf.readline()
        atom_index, segment, resid, resname, atomname, atomtype, charge, atomicweight, junk = l.strip().split()
        resid = int(resid)
        charge = float(charge)

        if resname == "CYS" and atomname == "SG":
            cys_sg_list.append(atom_index)
        
        charges.append(charge)
        atom_types.append(atomtype)

        # Are we in a new residue? If so, increment the residue count
        if last_resid != resid or last_segment != segment:
            num_residues += 1
        last_resid, last_segment = resid, segment
        
        if found_end_of_solute is False and resname == "TIP3":
            found_end_of_solute = True
            num_solute_residues = num_residues - 1 # We subtract 1 because the residue count has incremented into the first non-solute atom
            print >>sys.stderr, "Looks like you have %d atoms and %d residues in the solute (guessed by taking the atoms before the first TIP3)." % (num_solute_atoms, num_solute_residues)
            
        if found_end_of_solute is False:
            num_solute_atoms += 1
            # Unlike the AMBER residue map, ours is zero-based
            residue_map.append(num_residues - 1)

    psf.readline() # Blank line
    print >>sys.stderr, "%s contains %d residues as far as I can tell" % (filename, num_residues)

    if found_end_of_solute is False:
        num_solute_residues = num_residues

    # XXX: We don't care about the bonded terms or whatever, for this, so we stop reading

    # OK. Now we're done reading the PSF file.    
    return atom_types, charges


def load_prm(filename):
    prm = open(filename, 'r')
    # We need to extract the LJ terms from the prm file. Happily, we don't
    # need anything else from it (since we are only doing nonbonded energies).
    l = ""
    params = {}
    while l[0:4] != "NONB":
        l = prm.readline()
        if l == '':
            print >>sys.stderr, 'load_prm: Seems to be no NONB block in %s. That is probably OK.' % filename
            return params
        l = l.strip()
    if l[-1] == '-': l = prm.readline().strip()  # Eat extra nonsense

    nuke_comment_re = re.compile('!.*$')  # Regexp to get rid of comments
    # We assume the HBOND line is the end of the nonbonded parameters

    while True:
        l = prm.readline()
        if l == "": break  # EOF
        # Skip continuation line
        if l.strip().endswith('-'):
            prm.readline()
            continue
        l = nuke_comment_re.sub('', l).strip()
        if l == "": continue
        # Less than 4 tokens means this is probably the end of the section
        if len(l.split()) < 4: break
        if l[0:4] == "HBON" or l[0:4] == "NBFI": break
        # If the line doesn't start with a space, it's a new section, and we should stop
        # TODO: What's NBFIX?
        atomtype = ""
        epsilon, rmin2, eps14, rmin214 = None, None, None, None
        try:
            (atomtype, dummy, epsilon, rmin2, dummy, eps14, rmin214) = l.split()
            eps14 = float(eps14.strip())
            rmin214 = float(rmin214.strip())
        except ValueError:
            (atomtype, dummy, epsilon, rmin2) = l.split()
        epsilon = float(epsilon.strip())
        rmin2 = float(rmin2.strip())
        # Put the parameters in a dict for easy consumption
        this_param = {}
        this_param['epsilon'] = epsilon
        this_param['rmin2'] = rmin2
        this_param['eps14'] = eps14
        this_param['rmin214'] = rmin214
        params[atomtype.strip()] = this_param

    return params


# FFTK moves the water closer by -0.2 A (default) from the QM coordinates
# for a "better approximation" during the MM calculations. So that's what
# we do here.
def move_water(coords, offset=-0.2):
    def distance_between_two_atoms(coords, i, j):
        _, xi, yi, zi = coords[i]
        _, xj, yj, zj = coords[j]
        dist = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
        return dist ** 0.5

    # Find the H in the water that is the closest to any part of the ligand
    best_dist = 9999999.0
    best_water_atom = None
    best_ligand_atom = None
    for water_atom in range(len(coords) - 3, len(coords)):
        atomic_num, x, y, z = coords[water_atom]
        if atomic_num == 1: # Is this a hydrogen?
            # If so, calculate the distance between this hydrogen and all atoms of the ligand,
            # keeping a note of the closest ligand atom
            for ligand_atom in range(len(coords) - 3):
                dist = distance_between_two_atoms(coords, water_atom, ligand_atom)
                if dist < best_dist:
                    best_dist = dist
                    best_ligand_atom = ligand_atom
                    best_water_atom = water_atom

    # print 'move_water: Closest distance is %.3f, between atoms %d and %d' % (best_dist, best_ligand_atom, best_water_atom)

    # Determine a vector of norm == offset and pointing toward the closest ligand atom
    _, lx, ly, lz = coords[best_ligand_atom]
    _, wx, wy, wz = coords[best_water_atom]
    dx, dy, dz = lx - wx, ly - wy, lz - wz
    # Normalize this vector to turn it into a unit vector
    dx, dy, dz = dx/best_dist, dy/best_dist, dz/best_dist
    dx, dy, dz = dx * offset, dy * offset, dz * offset

    shifted_coords = coords
    for water_atom in range(len(coords) - 3, len(coords)):
        atomic_num, x, y, z = coords[water_atom]
        x, y, z = x + dx, y + dy, z + dz
        shifted_coords[water_atom] = atomic_num, x, y, z

    # Return a coords array where the water only has been displaced by that vector
    return shifted_coords


def calc_mm_interaction_energy(coords, water_shift, ligand_psf, ff_prm):
    # TODO:
    # Assert that the last three atoms' atomic numbers are 1, 8, 1
    # assert(...)

    # Hardcoded MM parameters for water - charges and L-J
    # Dicts are indexed by atom type
    tip3p_charges = {1: 0.417, 8: -0.834}
    tip3p_eps = {1: -0.046, 8: -0.1521}
    tip3p_rmin = {1: 0.2245, 8: 1.7682}
    # qH, qO = 0.417, -0.834
    # epsH, epsO = -0.046, -0.1521
    # rminH, rminO = 0.2245, 1.7682

    coords = move_water(coords, water_shift)

    ligand_atom_types, ligand_charges = ligand_psf

    total_electro, total_lj = 0.0, 0.0
    # Again we assume that the last three atoms are a water molecule
    for i in range(len(coords) - 3):
        atomtype = ligand_atom_types[i]
        charge = ligand_charges[i]
        atomic_num, x, y, z = coords[i]
        # Iterate over each water atom
        electro, lj = 0.0, 0.0
        for j in range(len(coords) - 3, len(coords)):
            tip3p_atomic_num, tip3p_x, tip3p_y, tip3p_z = coords[j]
            # Calculate euclidean distance between this water particle and the ligand atom
            dist = (x - tip3p_x)**2 + (y - tip3p_y)**2 + (z - tip3p_z)**2
            dist = dist**0.5

            # Electrostatic energy
            electro += COULOMB * ((charge * tip3p_charges[tip3p_atomic_num]) / dist)
            # Lennard-Jones energy
            lj6 = ((tip3p_rmin[tip3p_atomic_num] + ff_prm[atomtype]['rmin2']) / dist)**6
            lj += (tip3p_eps[tip3p_atomic_num] * ff_prm[atomtype]['epsilon'])**0.5 * (lj6**2 - 2*lj6)

        total_electro += electro
        total_lj += lj

    # print 'MM total interaction energy of ligand with TIP3P = %.3f' % (total_electro + total_lj)
    return total_electro, total_lj



def main():
    ap = argparse.ArgumentParser(description='Compare QM vs MM water interaction energies')
    # Default scale value taken from isoflurane paper and its reasoning about discrepancy between
    # QM B3LYP water dimer interaction energy and MM TIP3P water dimer interaction energy
    ap.add_argument('--scale', type=float, default=-6.95/-4.84, help='Multiply QM energies by this number before printing them')
    ap.add_argument('--shift', type=float, default=0.0, help='Change the distance of the TIP3P water from ligand by this many Angstroms before calculating MM energy')
    ap.add_argument('system_yaml')
    ap.add_argument('waterint_logs', nargs='+')
    args = ap.parse_args()
    system = yaml.load(open(args.system_yaml))

    ff_prm = {}
    for filename in system['prms']:
        print >>sys.stderr, 'Loading force field parameter file %s...' % filename
        ff_prm.update(load_prm(filename))

    print >>sys.stderr, 'After loading all force field parameter files, I know about %d atomtypes.' % len(ff_prm.keys())

    print >>sys.stderr, 'Loading ligand PSF %s...' % system['psf']
    ligand_psf = load_psf(system['psf'])

    # We use single-point energies for the compound in question as well as water
    # as a baseline. For some reason FFTK only uses Hartree-Fock.
    # Will uncomment two lines below when we have some B3LYP energies
    sp_energies, _, _ = parse_waterint_log(system['sp_ligand_log'])
    sp_wat_energies, _, _ = parse_waterint_log(system['sp_wat_log'])

    # Extract the final QM energies from the water interaction calculations.
    print >>sys.stderr, 'Scale factor for QM energies: %.4f.' % args.scale
    print >>sys.stderr, 'Moving MM water by this much: %.2f Angstroms (negative means closer to ligand).' % args.shift
    print >>sys.stderr, ''
    print >>sys.stderr, 'Single-point ligand energy: %.4f kcal/mol' % sp_energies[-1]
    print >>sys.stderr, 'Single-point water energy: %.4f kcal/mol' % sp_wat_energies[-1]
    print >>sys.stderr, 'Interaction energies - single TIP3P vs ligand:'

    out = None

    for filename in args.waterint_logs:
        scf_energies, qm_coords, qm_params = parse_waterint_log(filename)

        for frame_i in range(len(qm_coords)):
            # Only show optimized energies
            if qm_params[frame_i] == {}:
                continue
            # Don't fire up the CSV writer until we need it, because we don't know all the field names until now
            if out is None:
                fieldnames = ['filename', 'qm_energy', 'mm_electro', 'mm_lj', 'mm_energy']
                fieldnames.extend(qm_params[frame_i].keys())
                out = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
                out.writeheader()

            data = qm_params[frame_i]
            data['filename'] = filename
            data['qm_energy'] = args.scale * (scf_energies[frame_i] - sp_energies[-1] - sp_wat_energies[-1])
            data['mm_electro'], data['mm_lj'] = calc_mm_interaction_energy(qm_coords[frame_i], args.shift, ligand_psf, ff_prm)
            data['mm_energy'] = data['mm_electro'] + data['mm_lj']
            out.writerow(data)


if __name__ == '__main__':
    sys.exit(int(main() or 0))
