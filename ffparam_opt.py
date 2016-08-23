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
from MDAnalysis.core.topologyobjects import Dihedral

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
                atom.occupancy = 1.0
                atom.tempfactor = 1.0
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
        # if tokens[0] == 'ETITLE:': print(line)
        if tokens[0] == 'ENERGY:': energy = float(tokens[energy_index])
        # if tokens[0] == 'ENERGY:' and tokens[1] == '1000': print(line)

    return energy


def calc_one_mm_energy(data):
    # Load template PDB
    pdb = PDB(data['pdb'])
    outputname = 'Dihedral_' + '_'.join([str(x) for x in data['indices']]) + '_' + str(data['dihedral_angle'])
    
    # Replace coordinates with the supplied ones and write out a new (temporary) PDB
    coords = data['coords']
    for i in range(len(pdb.atoms)):
        pdb.atoms[i].x = coords[i][0]
        pdb.atoms[i].y = coords[i][1]
        pdb.atoms[i].z = coords[i][2]
    pdb_fname = '%s.initial.pdb' % outputname
    pdb_fp = open(pdb_fname, 'w')   
    pdb.write(pdb_fp)
    pdb_fp.close()
    
    # Extrabonds file is used to zero out the dihedrals we are fitting
    # To do this, we need to associate atom types to atom indices, perhaps with the help of MDAnalysis
    u = mda.Universe(data['psf'], pdb_fname)
    dihedral_atomtypes = []
    for i in data['indices']: dihedral_atomtypes.append(u.atoms[i].type)
    data['dihedral_atomtypes'] = dihedral_atomtypes
    # DEBUG
    pdb.write(open('%s.initial.pdb' % outputname, 'w'))

    # Construct parameters block
    prm_string = ''
    for prm in data['prms']:
        prm_string += 'parameters %s/%s\n' % (getcwd(), prm)

    # Write the dihedrals to this file
    zero_prm_fname = '%s.zero.prm' % outputname
    zero_prm_fp = open(zero_prm_fname, 'w')
    extrabonds_fname = '%s.extrabonds.txt' % outputname
    extrabonds_fp = open(extrabonds_fname, 'w')
    
    zero_prm_fp.write('BONDS\nANGLES\nDIHEDRAL\n')
    for indices in data['all_dihedrals_to_be_fit']:
        names = []
        for i in data['indices']: names.append(u.atoms[i].type)
        #for multiplicity in (1, 2, 3, 4, 6):
        #    zero_prm_fp.write('%s 0.0 %d 0.0\n' % (' '.join(names), multiplicity))
        # Measure the dihedral angle for this dihedral, and put an entry in
        # extrabonds.txt that fixes that dihedral in place using a very high force constant
        these_atoms = [u.atoms[i] for i in indices]
        dihedral_angle = Dihedral(these_atoms).value()
        extrabonds_fp.write('dihedral %d %d %d %d 10000. %f\n' % (
                    indices[0], indices[1], indices[2], indices[3],
                    dihedral_angle))
    zero_prm_fp.write('\nEND\n')
    zero_prm_fp.close()
    extrabonds_fp.close()
    prm_string += 'parameters %s\n' % (zero_prm_fname)

    extrabonds_string = ''
    # Write stuff to extrabonds if you want
    extrabonds_string = 'extraBonds yes\nextraBondsFile %s' % extrabonds_fname
    
    # First, we run NAMD to minimize the structure with the dihedrals in question fixed using a large
    # force constant. Then, we run it *again*, with slightly different parameters, to actually calcualte
    # the molecular mechanics energy. On this second run, we don't use an extrabonds.txt because
    # we just want a single point energy. But we do keep the dihedral parameters that we are fitting
    # zeroed out in this energy calculation, because we would want to calculate the dihedral contributions
    # ourselves, so the user knows how much energy is being contributed by manipulating the dihedral parameters.
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
binaryrestart no
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
    namdconf_fname = '%s.namd.min.conf' % outputname
    namdconf_fp = open(namdconf_fname, 'w')
    namdconf_fp.write(namd_conf)
    namdconf_fp.close()

    namd_wd = dirname(abspath(namdconf_fname))   
    p = subprocess.Popen([data['namd'], namdconf_fname], stdout=subprocess.PIPE, cwd=getcwd())
    (namdlog, namdstderr) = p.communicate()

    # Now run NAMD a second time to calculate a single point energy
    # using the minimized coordinates from the previous NAMD run.
    namd_conf = """structure %(psf)s
coordinates %(pdb_fname)s
paraTypeCharmm on
%(prm_string)s
numsteps 1
exclude scaled1-4
outputName %(outputname)s
temperature 0
COMmotion yes
dielectric 1.0
cutoff 1000.0
switchdist 10.0
pairInteraction		 on
pairInteractionGroup1 1
pairInteractionFile   %(pdb_fname)s
pairInteractionSelf		 on
run 0
""" % {'psf': "%s/%s" % (getcwd(), data['psf']),
            'prm_string': prm_string,
            'pdb_fname': '%s.restart.coor' % outputname,
            'outputname': outputname}
    
    namdconf_fname = '%s.namd.energy.conf' % outputname
    namdconf_fp = open(namdconf_fname, 'w')
    namdconf_fp.write(namd_conf)
    namdconf_fp.close()
    p = subprocess.Popen([data['namd'], namdconf_fname], stdout=subprocess.PIPE, cwd=getcwd())
    (namdlog, namdstderr) = p.communicate()

    # DEBUG
    open(outputname + '.namd2.log', 'w').write(namdlog)
    
    # TODO: Harvest the dihedral angles we are trying to fit
    data['mm_energy'] = get_energy_from_namd_log(namdlog, 'TOTAL')
    print('%s at %.1f deg: MM=%.2f kcal/mol, QM=%.2f kcal/mol' % (' '.join(dihedral_atomtypes),
                data['dihedral_angle'], 
                data['mm_energy'], 
                data['qm_energy']))


    # Delete a file without caring whether it succeeds, such as if the file didn't exist in the first place
    def unlink_dontcare(fname):
        try:
            # unlink(fname)
            pass # DEBUG
        except:
            pass # YOLO

    # Remove the temporary files we made
    for fname in (pdb_fname, extrabonds_fname, namdconf_fname):
        unlink_dontcare(fname)

    # Remove NAMD output that we don't care about
    # The try/except blocks don't do anything because we don't care about failure
    for suffix in ('.coor', '.vel', '.xsc'):
        for suffix2 in ('', '.old', '.BAK'):
            unlink_dontcare('%s/%s%s%s' % (namd_wd, outputname, suffix, suffix2))
            unlink_dontcare('%s/%s.restart%s%s' % (namd_wd, outputname, suffix, suffix2))

    return data

def calc_mm_energy(psf, pdb, prms, dihedral_results, namd='namd2'):
    # For each entry in dihedral_results, calculate a relaxed MM energy.
    # This involves fixing all the dihedrals that we are trying to fit, and letting
    # NAMD minimize the rest of the structure. We use the input coordinates from
    # the Gaussian log files.
    # Use the CHARMM force field parameter files in prms.
    # We can also run a bunch of NAMDs simultaneously to make this faster.

    all_dihedrals_to_be_fit = []
    for d in dihedral_results:
        if d['indices'] not in all_dihedrals_to_be_fit:
            all_dihedrals_to_be_fit.append(d['indices'])

    print(all_dihedrals_to_be_fit)

    data = []
    for i in range(len(dihedral_results)):
        d = {
            'index': i,
            'psf': psf,
            'pdb': pdb,
            'prms': prms,
            'namd': namd,
            'all_dihedrals_to_be_fit': all_dihedrals_to_be_fit
        }
        d.update(dihedral_results[i])
        data.append(d)

    p = Pool(12)
    data = p.map(calc_one_mm_energy, data)
    #for i in range(len(data)):
    #    data[i] = calc_one_mm_energy(data[i])

    def shift_zeropoint(data, key):
        smallest = data[0][key]
        for d in data: smallest = d[key] if d[key] < smallest else smallest
        for d in data: d[key] -= smallest


    shift_zeropoint(data, 'mm_energy')
    shift_zeropoint(data, 'qm_energy')

    return data

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
            for j in range(1, 5): indices.append(int(matched.group(j)) - 1)

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
                                     'qm_energy': current_energy,
                                     'coords': coords})

        i += 1 # Advance line counter
    return dihedral_results


def results_to_html(data, fname='results.html'):
    f = open(fname, 'w')

    qm_energy, mm_energy, x_vals = [], [], []
    for x in range(len(data)):
        x_vals.append(x)
        qm_energy.append(data[x]['qm_energy'])
        mm_energy.append(data[x]['mm_energy'])

    qm_energy_str = 'x: %s, y: %s' % (x_vals, qm_energy)
    mm_energy_str = 'x: %s, y: %s' % (x_vals, mm_energy)
    
    frame_info_str = """<table class='table table-condensed table-hover table-striped'>
    <tr><th>Frame</th><th>Dihedral being scanned</th><th>Atom indices of that dihedral</th><th>Angle (degrees)</th></tr>"""
    for i in range(len(data)):
        d = data[i]
        frame_info_str += '<tr><td>%d</td><td>%s</td><td>%s</td><td>%.2f</td></tr>' % (
            i, ' '.join(d['dihedral_atomtypes']), d['indices'], d['dihedral_angle'])
    frame_info_str += "</table>"

    f.write(
"""<html><head>
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css">
<script src="https://code.jquery.com/jquery-3.1.0.slim.min.js"></script>
<script src="https://cdn.plot.ly/plotly-1.2.0.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>
<title>Dihedral energy surface</title>
</head>
<body>
<div class="container">
<div id="energy_plot" style="width: 100%%; height: 500px;"></div>
<div id="info_table">%s</div>
</div>
<script>
var qm_energy = { %s, name: 'QM' };
var mm_energy = { %s, name: 'MM (with dihedrals)' };
var data = [qm_energy, mm_energy];

$(document).ready(function() {
    var layout = {showlegend: true, legend: {"orientation": "h"},
        xaxis: {
            title: 'Frame'
        },
        yaxis: {
            title: 'Energy (kcal/mol)'
        }
    };
    Plotly.newPlot('energy_plot', data, layout);
});

</script>
</body>
</html>
""" % (frame_info_str, qm_energy_str, mm_energy_str))

    f.close()


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Something with force fields or whatever')
    ap.add_argument('system_yaml')
    ap.add_argument('dihedral_scan_logs', nargs='+')
    args = ap.parse_args()
    system = yaml.load(open(args.system_yaml))
    # We want to know what the QM energies that we are targeting are.
    # Hence we read in the hilariously unstructured output of Gaussian.
    dihedral_results = []
    print('Loading %d Gaussian dihedral scan logs.' % len(args.dihedral_scan_logs))
    for fname in args.dihedral_scan_logs:
        dihedral_results.extend(parse_gaussian_dihedral_scan_log(fname))

    # Now we calculate the relaxed MM energy for each of those conformations.
    data = calc_mm_energy(system['psf'], system['pdb'], system['prms'], dihedral_results)

    # Finally we should be able to make a plot of these energies, and make a table that
    # associates these energies to the actual dihedral angle we varied. Amazing!
    results_to_html(data)

    # End goal here is to be able to dynamically edit the PRM file we are creating
    # for this ligand, and watch the MM energy surface change in relation to the QM
    # energy surface. Should in the end speed up this whole process.
