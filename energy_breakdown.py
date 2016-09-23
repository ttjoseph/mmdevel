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
        for section in ('ATOM', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'NONB'):
            if tokens[0].startswith(section):
                current_section = section
                skip_line = True
                break

        if skip_line: continue

        if current_section is 'BOND':
            key = '%s %s' % (tokens[0], tokens[1])
            ffp.bonds[key] = {
                'force_constant': float(tokens[2]),
                'equilibrium_distance': float(tokens[3])
            }
        elif current_section is 'ANGL':
            # TODO: Urey-Bradley terms
            key = '%s %s %s' % (tokens[0], tokens[1], tokens[2])
            ffp.angles[key] = {
                'force_constant': float(tokens[3]),
                'equilibrium_angle': float(tokens[4]) * pi / 180.0
            }
        elif current_section is 'DIHE':
            key = '%s %s %s %s' % (tokens[0], tokens[1], tokens[2], tokens[3])
            ffp.dihedrals[key] = {
                'force_constant': float(tokens[4]),
                'multiplicity': float(tokens[5]),
                'delta': float(tokens[6])
            }
        elif current_section is 'IMPR':
            key = '%s %s %s %s' % (tokens[0], tokens[1], tokens[2], tokens[3])
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

    # Evaluate each bond term's contribution
    # V(bond) = Kb(b - b0)**2
    for group in u.bonds:
        key = ' '.join(group.type)
        for ff in force_field:
            if key in ff.bonds:
                term = ff.bonds[key]
                energy = term['force_constant'] * (group.length() - term['equilibrium_distance'])**2
                print('[%s] Bond: %s - %s: %.2f' % (ff.name, get_group_label(group), key, energy))
                if energy > 200:
                    print('length: %f, equil_length: %f' % (group.length(), term['equilibrium_distance']))
 
    # Evaluate each force field term's contribution to the angle energy
    # V(angle) = Ktheta(Theta - Theta0)**2
    # TODO: Calculate Urey-Bradley contributions
    # V(Urey-Bradley) = Kub(S - S0)**2
    for group in u.angles:
        key = ' '.join(group.type)
        for ff in force_field:
            if key in ff.angles:
                term = ff.angles[key]
                energy = ((group.angle() * pi / 180.0) - term['equilibrium_angle'])**2
                energy *= term['force_constant']
                print('[%s] Angle: %s - %s: %.2f' % (ff.name, get_group_label(group), key, energy))

    # And the same deal for dihedrals
    for group in u.dihedrals:
        key = ' '.join(group.type)
        for ff in force_field:
            if key in ff.dihedrals:
                term = ff.dihedrals[key]
                # force_constant, multiplicity, delta
                # V(dihedral) = Kchi(1 + cos(n(chi) - delta))
                energy = term['force_constant'] * (1 + cos(term['multiplicity'] * (
                    group.dihedral() * pi / 180.0) - term['delta']))
                print('[%s] Dihedral: %s - %s: %.2f' % (ff.name, get_group_label(group), key, energy))


if __name__ == '__main__':
    main()