#!/usr/bin/env python
import sys
import argparse
import re
import yaml
from psf_prm_to_bin import get_psf_block_header, read_int_block

# Constants lifted from FFTK
HARTREES_PER_KCAL = 1.041308e-21
THINGS_PER_MOLE = 6.02214e23
COULOMB = 332.0636
QM_ENERGY_SCALE = 1.16 # From FFTK, and currently a mystery

def parse_waterint_log(filename):
    """Extracts SCF energies from a water interaction calculation done by Gaussian."""
    lines = open(filename).readlines()
    energies, coords = [], []
    i = 0
    while i < len(lines):
        if re.search(r'Error termination', lines[i]):
            sys.stderr.write('Careful, %s has a Gaussian error termination.' % filename)
            return energies
        
        if re.search(r'Normal termination of Gaussian', lines[i]):
            # All is well
            # sys.stderr.write('Gaussian did not complain in %s. Nice!\n' % filename)
            break
        
        # Here, FFTK parses a relaxed potential scan...do we need to?
        # XXX: For now we'll not worry about it

        if re.search(r'Standard orientation:', lines[i]):
            # Save the input coordinates for this piece of the dihedral scan
            # Eat the header
            i += 5
            # Read the coordinates
            current_frame_coords = []
            while lines[i][1] != '-':
                (idx, atomic_num, atomic_type, x, y, z) = lines[i].split()
                x, y, z = float(x), float(y), float(z)
                current_frame_coords.append((int(atomic_num), x, y, z))
                i += 1
            coords.append(current_frame_coords)

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
            energies.append(scf)
            # print '%s: Energy: %f (relative: %f)' % (filename, scf, scf - energies[0])

        i += 1
    
    return energies, coords


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
    print "Number of atoms: %d" % num_atoms
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
            print "Looks like you have %d atoms and %d residues in the solute (guessed by taking the atoms before the first TIP3)." % (num_solute_atoms, num_solute_residues)
            
        if found_end_of_solute is False:
            num_solute_atoms += 1
            # Unlike the AMBER residue map, ours is zero-based
            residue_map.append(num_residues - 1)

    psf.readline() # Blank line
    print "%s contains %d residues as far as I can tell" % (filename, num_residues)

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
            print 'Seems to be no NONB block in %s. That is probably OK.' % filename
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


def calc_mm_interaction_energy(coords, ligand_psf, ff_prm):
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
    ap.add_argument('system_yaml')
    ap.add_argument('waterint_logs', nargs='+')
    args = ap.parse_args()
    system = yaml.load(open(args.system_yaml))

    ff_prm = {}
    for filename in system['prms']:
        print 'Loading force field parameter file %s...' % filename
        ff_prm.update(load_prm(filename))

    print 'I know about %d atomtypes now.' % len(ff_prm.keys())

    print 'Loading ligand PSF %s...' % system['psf']
    ligand_psf = load_psf(system['psf'])

    # We use single-point energies for the compound in question as well as water
    # as a baseline. For some reason FFTK only uses Hartree-Fock.
    sp_hf_energies, _ = parse_waterint_log(system['sp_hf_log'])
    sp_mp2_energies, _ = parse_waterint_log(system['sp_mp2_log'])
    sp_wat_energies, _ = parse_waterint_log(system['sp_wat_log'])

    # Extract the final QM energies from the water interaction calculations.
    print '\nInteraction energies - single TIP3P vs ligand:'
    energies = {}
    for filename in args.waterint_logs:
        scf_energies, qm_coords = parse_waterint_log(filename)
        energies[filename] = {}
        energies[filename]['qm'] = scf_energies[-1] - sp_hf_energies[-1] - sp_wat_energies[-1]
        energies[filename]['qm'] *= QM_ENERGY_SCALE
        energies[filename]['electro'], energies[filename]['lj'] =  calc_mm_interaction_energy(qm_coords[-1], ligand_psf, ff_prm)

        print '%s: QM = %.3f, MMelec = %.3f, MMlj = %.3f, MMtot = %.3f' % (filename, energies[filename]['qm'],
                                            energies[filename]['electro'],
                                            energies[filename]['lj'],
                                            energies[filename]['electro'] + energies[filename]['lj'])
    
    # At this point we would want to compare each QM energy to the corresponding MM energy.
    # We now have the interaction energy between the water molecule and the compound in question.
    # At this point, nonbonded terms don't matter because they have in principle been canceled out.
    # So, just the regular old MM electrostatics should be good enough?

if __name__ == '__main__':
    sys.exit(int(main() or 0))
