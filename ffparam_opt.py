#!/usr/bin/env python
#
# Helper for optimizing force field parameters using NAMD
import argparse
import re
from multiprocessing import Pool
import tempfile
import sys
from os import fdopen, unlink, getcwd
from os.path import abspath, dirname
import subprocess
import yaml
import MDAnalysis as mda

class AtomRecord:
    '''Represents a single atom.'''

    def __init__(self, line):
        self.atomid = int(line[6:11])
        self.atomname = line[12:16].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21]
        if self.chain == " ":
            self.chain = "X"
        # Solvent residue ids are frequently mangled
        # so just set the mangled ones to 0
        try:
            self.resid = int(line[22:26])
        except:
            self.resid = 0
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = float(line[54:60])
        self.tempfactor = float(line[60:66])


class PDB:
    solvent_residues = ["WAT", "T3P", "HOH"]

    def __init__(self, filename):
        self.atoms = []
        self.load_from_pdb(filename)

    def load_from_pdb(self, filename):
        pdb = open(filename)
        for line in pdb:
            if line[0:6] == "ATOM  ":
                atom = AtomRecord(line)
                self.atoms.append(atom)
        # print >> sys.stderr, "Read %d atoms from %s." % (len(self.atoms), filename)

    def nuke_solvent(self):
        print >> sys.stderr, "Stripping solvent residues."
        self.atoms = [atom for atom in self.atoms if atom.resname not in self.solvent_residues]

    def write(self, fp):
        """Writes this PDB back out, in PDB format."""

        for atom in self.atoms:
            fp.write("ATOM  %5d %4s %3s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % \
                     (atom.atomid, atom.atomname, atom.resname, atom.resid,
                      atom.x, atom.y, atom.z,
                      atom.occupancy, atom.tempfactor))

def get_energy_from_namd_log(namdlog, energy_type='TOTAL'):
    """Parse a string containing a NAMD log from a single point relaxed MD simulation and
    extract the total energy"""
    energy_index = None
    energy = None
    for line in namdlog.split('\n'):
        tokens = line.split()
        if len(tokens) == 0: continue
        if tokens[0] == 'ETITLE:': energy_index = tokens.index(energy_type)   
        if tokens[0] == 'ENERGY:': energy = float(tokens[energy_index])

    return energy


def calc_one_mm_energy(data):
    # Load template PDB
    pdb = PDB(data['pdb'])

    # Replace coordinates with the supplied ones and write out a new (temporary) PDB
    coords = data['coords']
    for i in range(len(pdb.atoms)):
        pdb.atoms[i].x = coords[i][0]
        pdb.atoms[i].y = coords[i][1]
        pdb.atoms[i].z = coords[i][2]
    (pdb_fp, pdb_fname) = tempfile.mkstemp(suffix='.pdb')
    pdb_fp = fdopen(pdb_fp, 'w')
    pdb.write(pdb_fp)
    pdb_fp.flush()

    # Construct parameters block
    prm_string = ''
    for prm in data['prms']:
        prm_string += 'parameters %s/%s\n' % (getcwd(), prm)

    # Extrabonds file is used to zero out the dihedrals we are fitting
    # TODO: To do this, we need to associate atom types to atom indices, perhaps with the help of MDAnalysis
    u = mda.Universe(data['psf'])
    dihedral_atomtypes = []
    for i in data['indices']:
        dihedral_atomtypes.append(u.atoms[i].type)
    # TODO: Write the dihedrals to this file
    extrabonds_string = ''
    (extrabonds_fp, extrabonds_fname) = tempfile.mkstemp('.extrabonds.txt')
    extrabonds_string = 'extraBonds yes\nextraBondsFile %s' % extrabonds_fname

    outputname = 'Dihedral_' + '_'.join(dihedral_atomtypes)

    # NAMD configuration for this run
    namd_conf = """structure %(psf)s
coordinates %(pdb_fname)s
paraTypeCharmm on
%(prm_string)s
temperature 310
exclude scaled1-4
1-4scaling 1.0
cutoff 1000.0
switching on
switchdist 1000.0
pairlistdist 1000.0
timestep 1.0
nonbondedFreq 2
fullElectFrequency 4
stepspercycle 20
outputName %(outputname)s
restartfreq 1000
%(extrabonds_string)s
minimize 1000
reinitvels 310
run 0
""" % {'psf': "%s/%s" % (getcwd(), data['psf']),
            'prm_string': prm_string,
            'extrabonds_string': extrabonds_string,
            'pdb_fname': pdb_fname,
            'outputname': outputname}

    # Actually run NAMD and harvest the energy results
    (namdconf_fp, namdconf_fname) = tempfile.mkstemp('.namd.conf')
    namdconf_fp = fdopen(namdconf_fp, 'w')
    namdconf_fp.write(namd_conf)
    namdconf_fp.flush()

    namd_wd = dirname(abspath(namdconf_fname))
    
    p = subprocess.Popen([data['namd'], namdconf_fname], stdout=subprocess.PIPE, cwd=getcwd())
    (namdlog, namdstderr) = p.communicate()

    # TODO: Harvest the dihedral angles we are trying to fit
    energy = get_energy_from_namd_log(namdlog)
    print('%s at %.1f deg: MM=%.2f kcal/mol, QM=%.2f kcal/mol' % (' '.join(dihedral_atomtypes), data['dihedral_angle'], energy, data['energy']))

    # Remove the temporary files we made
    for fname in (pdb_fname, extrabonds_fname, namdconf_fname):
        unlink(fname)

    # Remove NAMD output that we don't care about
    # The try/except blocks don't do anything because we don't care about failure
    for suffix in ('.coor', '.vel', '.xsc'):
        for suffix2 in ('', '.old', '.BAK'):
            try:
                unlink('%s/%s%s%s' % (namd_wd, outputname, suffix, suffix2))
            except:
                pass
            try:
                unlink('%s/%s.restart%s%s' % (namd_wd, outputname, suffix, suffix2))
            except:
                pass


def calc_mm_energy(psf, pdb, prms, dihedral_results, namd='namd2'):
    # For each entry in dihedral_results, calculate a relaxed MM energy.
    # This involves fixing all the dihedrals that we are trying to fit, and letting
    # NAMD minimize the rest of the structure. We use the input coordinates from
    # the Gaussian log files.
    # Use the CHARMM force field parameter files in prms.
    # We can also run a bunch of NAMDs simultaneously to make this faster.

    data = []
    for i in range(len(dihedral_results)):
        d = {
            'index': i,
            'psf': psf,
            'pdb': pdb,
            'prms': prms,
            'namd': namd
        }
        d.update(dihedral_results[i])
        data.append(d)

    p = Pool(8)
    p.map(calc_one_mm_energy, data)
# for d in data: calc_one_mm_energy(d)


def parse_gaussian_dihedral_scan_log(fname):
    """Parses Gaussian09 dihedral scan log file."""

    KCAL_MOL_PER_HARTREE = 627.5095

    def get_token(s, idx):
        tokens = s.split()
        if idx < len(tokens):
            return tokens[idx]
        else:
            return None

    with open(fname) as f:
        lines = f.readlines()

    i = 0
    dihedral_results = []
    while i < len(lines):
        # This signifies the start of a single point in the scan
        if re.search(r'Initial Parameters', lines[i]):
            # Chug along until we find the dihedral indices in question
            while get_token(lines[i], 4) != 'Scan': i += 1
            # ! D8    D(1,2,3,8)            177.9793         Scan
            # Use a horrifying regex to extract the dihedral indices
            indices_string = get_token(lines[i], 2)
            matched = re.match(r'.\((\d+),(\d+),(\d+),(\d+)\)', indices_string)
            indices = []
            for j in range(1, 5): indices.append(int(matched.group(j)))

        if re.search(r'Input orientation:', lines[i]):
            # Save the input coordinates for this piece of the dihedral scan
            # Eat the header
            i += 5
            # Read the coordinates
            coords = []
            while lines[i][1] != '-':
                (idx, atomic_num, atomic_type, x, y, z) = lines[i].split()
                x, y, z = float(x), float(y), float(z)
                coords.append((x, y, z))
                i += 1

        # This and the following blocks extract the current energy for this dihedral
        if re.search(r'SCF[ \t]*Done:', lines[i]):
            current_energy = float(get_token(lines[i], 4)) * KCAL_MOL_PER_HARTREE

        # We deliberately favor MP2 energy over the RHF energy extracted above
        if re.search(r'E2.*EUMP2', lines[i]):
            current_energy = float(get_token(lines[i], 5).replace('D', 'E')) * KCAL_MOL_PER_HARTREE

        if re.search(r'Optimization completed\.', lines[i]):
            # Eat lines until we find the line with our dihedral so we can extract its angle value
            while get_token(lines[i], 2) != indices_string: i += 1
            dihedral_angle = float(get_token(lines[i], 3))
            dihedral_results.append({'indices': indices,
                                     'dihedral_angle': dihedral_angle,
                                     'energy': current_energy,
                                     'coords': coords})

        i += 1 # Advance line counter
    return dihedral_results


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Something with force fields or whatever')
    ap.add_argument('system_yaml')
    ap.add_argument('dihedral_scan_logs', nargs='+')
    args = ap.parse_args()
    system = yaml.load(open(args.system_yaml))
    # We want to know what the QM energies that we are targeting are.
    # Hence we read in the hilariously unstructured output of Gaussian.
    dihedral_results = []
    for fname in args.dihedral_scan_logs:
        dihedral_results.extend(parse_gaussian_dihedral_scan_log(fname))

    # Now we calculate the relaxed MM energy for each of those conformations.
    calc_mm_energy(system['psf'], system['pdb'], system['prms'], dihedral_results)

    # Finally we should be able to make a plot of these energies, and make a table that
    # associates these energies to the actual dihedral angle we varied. Amazing!

    # End goal here is to be able to dynamically edit the PRM file we are creating
    # for this ligand, and watch the MM energy surface change in relation to the QM
    # energy surface. Should in the end speed up this whole process.
