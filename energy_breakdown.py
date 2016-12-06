#!/usr/bin/env python
#
# Breaks down the energy contributions of individual force field terms for a given structure.
# Useful for debugging messed up energies.
import re
import argparse
from math import pi, cos
import yaml
import MDAnalysis as mda

class ForceFieldParams(object):
    def __init__(self, name):
        self.name = name
        self.atomtypes = {}
        self.bonds = {}
        self.angles = {}
        self.dihedrals = {}
        self.impropers = {}
        self.cmap = {}

# Each force field parameter can have two names (e.g. A B C D and D C B A)
# so we need to care about that in order not to ignore existing parameters when
# checking for missing parameters 
def key_names(names):
    key1 = ' '.join(names)
    key2 = ' '.join(reversed(names))
    return key1, key2

def load_charmm_ff_params(fname):
    """Load a CHARMM .prm file containing force field parameters."""
    with open(fname) as f:
        lines = f.readlines()

    comment_stripper = re.compile(r'[!\*].*')
    ffp = ForceFieldParams(fname)

    current_section = None
    for i in range(len(lines)):
        # Ignore comments and blank lines
        line = comment_stripper.sub('', lines[i].strip())
        if line == '': continue

        tokens = line.split()
        skip_line = False
        for section in ('ATOM', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'NONB', 'CMAP'):
            if tokens[0].startswith(section):
                current_section = section
                skip_line = True
                break

        if skip_line: continue

        if current_section is 'BOND':
            key1, key2 = key_names((tokens[0], tokens[1]))
            ffp.bonds[key1] = ffp.bonds[key2] = {
                'force_constant': float(tokens[2]),
                'equilibrium_distance': float(tokens[3])
            }
        elif current_section is 'ANGL':
            # TODO: Urey-Bradley terms
            key1, key2 = key_names((tokens[0], tokens[1], tokens[2]))
            ffp.angles[key1] = ffp.angles[key2] = {
                'force_constant': float(tokens[3]),
                'equilibrium_angle': float(tokens[4]) * pi / 180.0
            }
        elif current_section is 'DIHE':
            key1, key2 = key_names((tokens[0], tokens[1], tokens[2], tokens[3]))
            ffp.dihedrals[key1] = ffp.dihedrals[key2] = {
                'force_constant': float(tokens[4]),
                'multiplicity': float(tokens[5]),
                'delta': float(tokens[6])
            }
        elif current_section is 'IMPR':
            key = key_names((tokens[0], tokens[1], tokens[2], tokens[3]))
        else:
            # Unknown line type
            continue
    return ffp

def get_group_label(group):
    """Make a human-readable label for an AtomGroup."""
    indices = [a.index for a in group.atoms]
    names = [a.name for a in group.atoms]
    label = []
    for i in range(len(indices)):
        label.append('%d/%s' % (indices[i], names[i]))
    return(' '.join(label))

def main():
    ap = argparse.ArgumentParser(description='Show contribution of individual force field terms to energies')
    ap.add_argument('system_yaml', help='YAML file describing the system. (PDB files listed in it are ignored.)')
    ap.add_argument('trajectory', help='Trajectory file. Single-structure PDB OK too. Used instead of whatever\'s in system.yaml')
    args = ap.parse_args()

    system = yaml.load(open(args.system_yaml))

    force_field = []

    for prm in system['prms']:
        force_field.append(load_charmm_ff_params(prm))

    u = mda.Universe(system['psf'], args.trajectory)

    unknown_param_keys = set()

    # Evaluate each bond term's contribution
    # V(bond) = Kb(b - b0)**2
    for group in u.bonds:
        key1, key2 = key_names(group.type)
        found_term = False
        for ff in force_field:
            if key1 in ff.bonds or key2 in ff.bonds:
                found_term = True
                term = ff.bonds[key1]
                energy = term['force_constant'] * (group.length() - term['equilibrium_distance'])**2
                print('[%s] Bond: %s - %s: %.2f' % (ff.name, get_group_label(group), key1, energy))
                if energy > 200:
                    print('length: %f, equil_length: %f' % (group.length(), term['equilibrium_distance']))
        if found_term is False:
            unknown_param_keys.add(key1)
 
    # Evaluate each force field term's contribution to the angle energy
    # V(angle) = Ktheta(Theta - Theta0)**2
    # TODO: Calculate Urey-Bradley contributions
    # V(Urey-Bradley) = Kub(S - S0)**2
    for group in u.angles:
        key1, key2 = key_names(group.type)
        found_term = False
        for ff in force_field:
            if key1 in ff.angles or key2 in ff.angles:
                found_term = True
                term = ff.angles[key1]
                energy = ((group.angle() * pi / 180.0) - term['equilibrium_angle'])**2
                energy *= term['force_constant']
                print('[%s] Angle: %s - %s: %.2f' % (ff.name, get_group_label(group), key1, energy))
        if found_term is False:
            unknown_param_keys.add(key1)

    # And the same deal for dihedrals
    for group in u.dihedrals:
        key1, key2 = key_names(group.type)
        found_term = False
        for ff in force_field:
            if key1 in ff.dihedrals or key2 in ff.dihedrals:
                found_term = True
                term = ff.dihedrals[key1]
                # force_constant, multiplicity, delta
                # V(dihedral) = Kchi(1 + cos(n(chi) - delta))
                energy = term['force_constant'] * (1 + cos(term['multiplicity'] * (
                    group.dihedral() * pi / 180.0) - term['delta']))
                print('[%s] Dihedral: %s - %s: %.2f' % (ff.name, get_group_label(group), key1, energy))
        if found_term is False:
            unknown_param_keys.add(key1)

    # TODO: eliminate redundant keys
    if len(unknown_param_keys) > 0:
        print('Could not find parameters for the following:')
        for x in sorted(unknown_param_keys):
            print(x)


if __name__ == '__main__':
    main()