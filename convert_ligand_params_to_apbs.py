#!/usr/bin/env python
#
# Want to include a custom ligand in an APBS calculation?
# Use this to convert the RTF and associated PRMs (in CHARMM force field format)
# to APBS force field format.
import sys
import argparse

# Parse the CHARMM force field parameters but just the NONB block, into a dict
def load_prm_nonb(prm_filename, atomtypes={}):
    in_nonb_block = False
    with open(prm_filename) as f:
        for line in f.readlines():
            l = line.strip()
            if l.startswith('NONB'):
                in_nonb_block = True
                continue
            elif l.startswith('!'):
                continue
            elif l.startswith('HBON'):
                break

            if in_nonb_block is True:
                try:
                    tokens = l.split()
                    atomname, eps_kcalmol, radius = tokens[0], tokens[2], tokens[3]
                    atomtypes[atomname] = {
                        'epsilon': float(eps_kcalmol) * 4.184,  # Convert to kJ/mol
                        'radius': float(radius)
                    }
                except Exception as e:
                    continue
    return atomtypes


# Loads an RTF
def load_rtf(rtf_filename):
    resname, charge, expected_charge = None, 0.0, None
    atoms = []
    with open(rtf_filename) as f:
        for line in f.readlines():
            l = line.strip()
            if l.startswith('ATOM'):
                tokens = l.split()
                atomname, atomtype, atomcharge = tokens[1], tokens[2], float(tokens[3])
                charge += atomcharge
                atoms.append({'name': atomname, 'type': atomtype, 'charge': atomcharge})
            elif l.startswith('RESI'):
                tokens = l.split()
                resname, expected_charge = tokens[1], float(tokens[2])

    print >> sys.stderr, 'Residue %s has expected charge %f and actual charge %f' % (resname, expected_charge, charge)
    return resname, atoms


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Convert CHARMM rtf/prm into APBS format')
    ap.add_argument('rtf', help='Residue topology file')
    ap.add_argument('prm', nargs='+', help='Force field parameter files')
    args = ap.parse_args()

    atomtypes = {}
    for prm in args.prm:
        atomtypes = load_prm_nonb(prm, atomtypes)

    print >> sys.stderr, 'Atom types parsed: %d' % len(atomtypes)

    resname, atoms = load_rtf(args.rtf)

    # Iterate through the ligand atoms and spit out the APBS force field block
    # resname atomname charge radius(A) epsilon(kJ/mol)
    print '# Residue %s from %s' % (resname, args.rtf)
    for atom in atoms:
        print '%s %s %f %f %f' % (resname,
                                  atom['name'].upper(), # PDB2PQR wants this uppercase for some reason
                                  atom['charge'],
                                  atomtypes[atom['type']]['radius'],
                                  atomtypes[atom['type']]['epsilon'])
