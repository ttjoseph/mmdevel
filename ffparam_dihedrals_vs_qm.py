#!/usr/bin/env python3
#
# Helper for optimizing force field parameters using NAMD
import argparse
import re
from multiprocessing import Pool, cpu_count
import tempfile
import sys
from os import fdopen, unlink, getcwd
from os.path import abspath, dirname, basename, splitext
import subprocess
from base64 import b64encode
import math
import yaml
import MDAnalysis as mda
from MDAnalysis.core.topologyobjects import Dihedral
from scipy.interpolate import griddata

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
        print("Stripping solvent residues.", file=sys.stderr)
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
    for line in namdlog.splitlines():
        tokens = line.split()
        if len(tokens) == 0: continue
        if tokens[0] == 'ETITLE:': energy_index = tokens.index(energy_type)
        if tokens[0] == 'ENERGY:': energy = float(tokens[energy_index])

    return energy


def calc_one_mm_energy(data, ignore_our_dihedrals=False):
    # Load template PDB and its associated PSF
    pdb = PDB(data['pdb'])

    # Make a string representing this data point, for human consumption
    description = []
    for key, dihedral in list(data['dihedrals'].items()):
        description.append(key + '_' + str(dihedral['angle']))
    description = '_'.join(description)

    # Pick a unique temporary file name
    outputname = f"Dihedral_{description}_{data['index']}"

    # Replace coordinates with the supplied ones and write out a new (temporary) PDB
    coords = data['coords']
    if len(pdb.atoms) != len(coords):
        print('PDB has %d atoms but I found %d atoms in the Gaussian output' % (len(pdb.atoms), len(coords)), file=sys.stderr)
    for i in range(len(pdb.atoms)):
        pdb.atoms[i].x = coords[i][0]
        pdb.atoms[i].y = coords[i][1]
        pdb.atoms[i].z = coords[i][2]
    pdb_fname = f'{outputname}.initial.pdb'
    pdb_fp = open(pdb_fname, 'w')
    pdb.write(pdb_fp)
    pdb_fp.close()

    # Extrabonds file contains the dihedrals we are fitting with very high force constants.
    # This holds those dihedrals static and forces the MM relaxation to happen around them.
    # To do this, we need to associate atom types to atom indices with the help of MDAnalysis
    u = mda.Universe(data['psf'], pdb_fname)
    for key, dihedral in list(data['dihedrals'].items()):
        dihedral['atomtypes'] = [u.atoms[i].type for i in dihedral['indices']]
        dihedral['atomnames'] = [u.atoms[i].name for i in dihedral['indices']]

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
    # Iterate over dihedrals
    for indices in data['all_dihedrals_to_be_fit']:
        names = []
        for i in indices: names.append(u.atoms[i].type)
        # We may want to calculate the relaxed molecular mechanics energy without contributions
        # from the dihedrals we are trying to optimize. In this case we could calculate
        # those separately and add them to the total energies. However in this case one
        # might expect the relaxed MM conformations to be slightly different from the 'real'
        # relaxed MM conformations with all energy terms included.
        if ignore_our_dihedrals:
            for multiplicity in (1, 2, 3, 4, 6):
                zero_prm_fp.write('%s 0.0 %d 0.0\n' % (' '.join(names), multiplicity))
        # Measure the dihedral angle for this dihedral, and put an entry in
        # extrabonds.txt that fixes that dihedral in place using a very high force constant
        dihedral_angle = Dihedral(indices, u).value()
        extrabonds_fp.write('dihedral %d %d %d %d 10000. %f\n' % (
                    indices[0], indices[1], indices[2], indices[3],
                    dihedral_angle))
    zero_prm_fp.write('\nEND\n')
    zero_prm_fp.close()
    extrabonds_fp.close()
    prm_string += 'parameters %s\n' % (zero_prm_fname)
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
    p = subprocess.Popen([data['namd'], namdconf_fname], stdout=subprocess.PIPE, cwd=getcwd(), universal_newlines=True)
    (namdlog, namdstderr) = p.communicate()

    # NAMD shouldn't be complaining. If it is, tell the user
    if namdstderr is not None:
        print(namdstderr)

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

    namdconf2_fname = '%s.namd.energy.conf' % outputname
    namdconf_fp = open(namdconf2_fname, 'w')
    namdconf_fp.write(namd_conf)
    namdconf_fp.close()
    p = subprocess.Popen([data['namd'], namdconf2_fname], stdout=subprocess.PIPE, cwd=getcwd(), universal_newlines=True)
    (namdlog, namdstderr) = p.communicate()
    namdlog = str(namdlog)

    # NAMD shouldn't be complaining. If it is, tell the user
    if namdstderr is not None:
        print(namdstderr)

    # Get individual energy types
    for energy_type in ('BOND', 'ANGLE', 'DIHED', 'IMPRP', 'ELECT', 'VDW'):
        data[energy_type] = get_energy_from_namd_log(namdlog, energy_type)
    data['mm_energy'] = get_energy_from_namd_log(namdlog, 'TOTAL')
    print(('%s: MM=%.2f kcal/mol, QM=%.2f kcal/mol' % (outputname,
                data['mm_energy'],
                data['qm_energy'])))


    # Delete a file without caring whether it succeeds, such as if the file didn't exist in the first place
    def unlink_dontcare(fname):
        try:
            unlink(fname)
        except:
            pass # YOLO

    # Remove the temporary files we made
    for fname in (pdb_fname, extrabonds_fname, namdconf_fname, namdconf2_fname):
        unlink_dontcare(fname)

    # Remove NAMD output that we don't care about
    for suffix in ('.coor', '.vel', '.xsc', '.zero.prm'):
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
    for result in dihedral_results:
        for dihedral_string in result['dihedrals']:
            dihedral = result['dihedrals'][dihedral_string]
            if dihedral['indices'] not in all_dihedrals_to_be_fit:
                all_dihedrals_to_be_fit.append(dihedral['indices'])

    print('Atom indices of dihedrals to be fit:')
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

    # Optimistically assume one core's worth of time will be spent in I/O starting up namd2 instances
    calc_one_mm_energy(data[0])
    p = Pool(cpu_count() + 1)
    data = p.map(calc_one_mm_energy, data)
    # For debugging only - stack trace in Python 2 is ruined by Pool
    # for i in range(len(data)):
    #     data[i] = calc_one_mm_energy(data[i])

    def shift_zeropoint(data, key):
        smallest = data[0][key]
        for d in data: smallest = d[key] if d[key] < smallest else smallest
        for d in data: d[key] -= smallest

    shift_zeropoint(data, 'mm_energy')
    shift_zeropoint(data, 'qm_energy')

    return data

def parse_gaussian_dihedral_scan_log(fname):
    """Parses Gaussian09 dihedral scan log file.
    
    Should now allow an arbitrary number of scanned (independently varying) dihedrals."""

    KCAL_MOL_PER_HARTREE = 627.5095

    def get_token(s, idx):
        tokens = s.split()
        if idx < len(tokens):
            return tokens[idx]
        else:
            return None

    # Read the entire file at once - these generally aren't that big
    with open(fname) as f:
        lines = f.readlines()

    i = 0
    dihedrals_scanned, dihedral_results = dict(), []
    dihedrals = None
    while i < len(lines):
        # This signifies the start of a single point in the scan
        if re.search(r'Initial Parameters', lines[i]):
            while lines[i].startswith(' GradGradGrad') is False:
                # Look for scanned dihedral and save it
                # ! D8    D(1,2,3,8)            177.9793         Scan
                if get_token(lines[i], 4) == 'Scan':
                    indices_string = get_token(lines[i], 2)
                    matched = re.match(r'.\((\d+),(\d+),(\d+),(\d+)\)', indices_string)
                    indices = []
                    for j in range(1, 5):
                        indices.append(int(matched.group(j)) - 1)
                    dihedrals_scanned[indices_string] = indices
                i += 1

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
            # We need to make a new dihedrals object for every point on the scan,
            # otherwise Python will keep reusing the same one and bork our results dict
            dihedrals = dict()
            for indices_string, indices in list(dihedrals_scanned.items()):
                dihedrals[indices_string] = {'indices': indices}

        # This and the following blocks extract the current energy for this dihedral set
        if re.search(r'SCF[ \t]*Done:', lines[i]):
            current_energy = float(get_token(lines[i], 4)) * KCAL_MOL_PER_HARTREE

        # We deliberately favor MP2 energy over the RHF or whatever energy extracted above
        if re.search(r'E2.*EUMP2', lines[i]):
            current_energy = float(get_token(lines[i], 5).replace('D', 'E')) * KCAL_MOL_PER_HARTREE

        # Save the result when we get to the final optimized coordinates, which we don't actually save
        # We also save each dihedral angle
        if re.search(r'Optimization completed\.', lines[i]):
            while lines[i].startswith(' GradGradGrad') is False and lines[i].startswith(' Iteration') is False:
                dihedral_string = get_token(lines[i], 2)
                if dihedral_string in dihedrals:
                    dihedrals[dihedral_string]['angle'] = float(get_token(lines[i], 3))
                i += 1
                if i >= len(lines):
                    print(f'We got problems with {fname}', file=sys.stderr)
            
            dihedral_results.append({'dihedrals': dihedrals,
                                     'gaussian_log_filename': fname,
                                     'qm_energy': current_energy,
                                     'coords': coords})

        i += 1 # Advance line counter

    return dihedral_results


def make_single_cmap_table(data, dihedral1, dihedral2, spacing=24, out=sys.stdout):
    """Generates a CMAP table for the given dihedral, given 2D dihedral scan data.
    
    We need to produce a 360x360 grid of energy terms that is [-180, 180) degrees for each dimension.
    The grid divides each dimension into equally-sized pieces, and we always start at -180 degrees.
    There is no CMAP "phase shift" to make fitting to the grid easier. And there are the same number
    of points in both dimensions.
    
    But we are only given an energy surface that covers part of the CMAP energy space, and it probably
    does not align to the preordained grid points. This is because the QM energy surface was
    probably generated with Gaussian09 Mod=OptRedundant, which does not seem to support pinning to
    specific dihedral angles.
    
    So, we construct the whole source energy grid, leaving the parts that were not sampled zero,
    then interpolate to map it onto the CMAP grid. The CMAP grid should probably have higher resolution
    than the QM energy surface, to minimize sampling error.
    """

    print(('! Here is the CMAP table for %s %s' % (dihedral1, dihedral2)))
    coords, values = [], []
    dihedral1_atomtypes, dihedral2_atomtypes = None, None
    for d in data:
        # If this is the dihedral pair in question, add it to our list of source data
        if dihedral1 in d['dihedrals'] and dihedral2 in d['dihedrals']:
            d1, d2 = d['dihedrals'][dihedral1], d['dihedrals'][dihedral2]
            dihedral1_atomtypes = d1['atomtypes']
            dihedral2_atomtypes = d2['atomtypes']
            # Use "periodic images" so griddata can interpolate at the edges of the
            # coordinate space
            for period1 in (-360, 0, 360):
                for period2 in (-360, 0, 360):
                    coords.append((d1['angle'] + period1, d2['angle'] + period2))
                    values.append(d['qm_energy'] - d['mm_energy'])
                    if (period1, period2) == (0, 0):
                        print(f"At {d1['angle']:.2f} {d2['angle']:.2f}, QM: {d['qm_energy']:.2f}, MM: {d['mm_energy']:.2f}, diff: {d['qm_energy'] - d['mm_energy']:.2f}",
                            file=sys.stderr)
    
    # This is probably a dumb way to do this
    xi = []
    for angle1 in range(-180, 180, int(math.floor(360/spacing))):
        for angle2 in range(-180, 180, int(math.floor(360/spacing))):
            xi.append((angle1, angle2))

    # We use scipy.interpolate.griddata. This takes coordinates (dihedral angles), values
    # (energy differences), and xi (CMAP grid coordinates), and returns interpolated values.
    cmap_raw = griddata(coords, values, xi, method='linear', fill_value=0.0)

    # Now to spit out the CMAP table so that NAMD/CHARMM will parse it
    out.write('CMAP\n')
    out.write('%s %s %d\n' % (' '.join(dihedral1_atomtypes), ' '.join(dihedral2_atomtypes), spacing))
    # Five values per line, why not
    for i in range(len(cmap_raw)):
        if i != 0 and i % 5 == 0: out.write('\n')
        out.write('   %.6f' % cmap_raw[i])
    out.write('\n\n')

def make_cmap_terms(data, fname=None):
    # Find all dihedral pairs and call make_single_cmap_table on each one.
    # Importantly, we tolerate the order of the dihedrals in each pair being reversed.
    # TODO: Does D(1,2,3,4) == D(4,3,2,1)?
    dihedral_pairs = set()
    for d in data:
        if len(d['dihedrals']) == 2:
            dihedral_pairs.add(' '.join(sorted(d['dihedrals'].keys())))

    out = open(fname, 'w') if fname is not None else sys.stdout
    
    for dihedral_pair in dihedral_pairs:
        dihedrals = dihedral_pair.split()
        make_single_cmap_table(data, dihedrals[0], dihedrals[1], out=out)
    

def generate_data_uri_download(fname):
    """Generate HTML that provides a download link for a file, where the file content
    is embedded as a data URI."""
    pretty_fname = basename(fname)
    with open(fname) as f:
        payload = b64encode(f.read().encode('utf8'))
        return '<a download="%s" href="data:text/plain;base64,%s">%s</a>' % (pretty_fname, payload, pretty_fname)


def results_to_html(data, fname='results.html'):
    f = open(fname, 'w')
    print(f'Writing HTML output to {fname}.')
    qm_energy, mm_energy, x_vals = [], [], []
    for x in range(len(data)):
        x_vals.append(data[x]['index'])
        qm_energy.append(data[x]['qm_energy'])
        mm_energy.append(data[x]['mm_energy'])

    output = {}
    output['qm_energy_str'] = 'x: %s, y: %s' % (x_vals, qm_energy)
    output['mm_energy_str'] = 'x: %s, y: %s' % (x_vals, mm_energy)

    # Make a pretty HTML table of what's shown in the plot
    frame_info_str = """<table class='table table-condensed table-hover table-striped'>
    <tr><th>Frame</th><th>Gaussian log</th><th>Atom types</th><th>Atom names</th><th>Atom indices</th>
    <th>Angle (degrees)</th></tr>"""
    for i in range(len(data)):
        indices, atomtypes, atomnames, angles = [], [], [], []
        for key, dihedral in list(data[i]['dihedrals'].items()):
            indices.append(dihedral['indices'])
            atomtypes.append(dihedral['atomtypes'])
            atomnames.append(dihedral['atomnames'])
            angles.append(dihedral['angle'])
        frame_info_str += '<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % (
            i, basename(data[i]['gaussian_log_filename']), 
            ' '.join(['-'.join(a) for a in atomtypes]),
            ' '.join(['-'.join(a) for a in atomnames]),
            ' '.join([str(a) for a in indices]),
            ' '.join([str(a) for a in angles]))
    frame_info_str += "</table>"
    output['frame_info_str'] = frame_info_str

    output['psf'] = data[0]['psf']
    output['pdb'] = data[0]['pdb']
    download_files = [output['psf'], output['pdb']]
    download_files.extend(data[0]['prms'])
    output['download_str'] = ' &sdot; '.join(map(generate_data_uri_download, download_files))
    f.write(
"""<html><head>
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css">
<script src="https://code.jquery.com/jquery-3.1.0.slim.min.js"></script>
<script src="https://cdn.plot.ly/plotly-1.2.0.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>
<title>Dihedral energy surface for %(psf)s</title>
<style>
#info_table {
  min-height: 400px;
  resize: both;
  overflow-x: hidden;
  overflow-y: auto;
  margin-bottom: 0.2em;
}

td {
    font-size: 75%%;
}

</style>
</head>
<body>
<div class="container-fluid">
<div id="energy_plot" style="width: 100%%; height: 500px;"></div>
<div id="info_table">%(frame_info_str)s</div>
Download files: %(download_str)s
</div>
<script>
var qm_energy = { %(qm_energy_str)s, name: 'QM' };
var mm_energy = { %(mm_energy_str)s, name: 'MM' };
var data = [qm_energy, mm_energy];

$(document).ready(function() {
    var layout = {
        showlegend: true,
        legend: {"orientation": "h"},
        title: 'Dihedral energy surface for %(psf)s',
        xaxis: { title: 'Frame' },
        yaxis: { title: 'Energy (kcal/mol)' }
    };
    Plotly.newPlot('energy_plot', data, layout);

    // Enforce minimum height of dihedral info table
    var newHeight = $(window).height() - $('.info_table').offset().top - 15;
    if(newHeight < 400) { newHeight = 400; }
    $('.info_table').height(newHeight);
});

</script>
</body>
</html>
""" % output)

    f.close()


def main():
    ap = argparse.ArgumentParser(description='Something with force fields or whatever')
    ap.add_argument('system_yaml')
    ap.add_argument('dihedral_scan_logs', nargs='+')
    args = ap.parse_args()
    system = yaml.load(open(args.system_yaml), Loader=yaml.SafeLoader)
    # We want to know what the QM energies that we are targeting are.
    # Hence we read in the hilariously unstructured output of Gaussian.
    dihedral_results = []
    print(f'Loading {len(args.dihedral_scan_logs)} Gaussian dihedral scan logs.')
    for fname in args.dihedral_scan_logs:
        dihedral_results.extend(parse_gaussian_dihedral_scan_log(fname))

    # Now we calculate the relaxed MM energy for each of those conformations.
    data = calc_mm_energy(system['psf'], system['pdb'], system['prms'], dihedral_results)

    make_cmap_terms(data, f'{system["psf"]}.cmap')
    print(f'Wrote CMAP terms to {system["psf"]}.cmap')

    # Finally we should be able to make a plot of these energies, and make a table that
    # associates these energies to the actual dihedral angle we varied. Amazing!
    results_to_html(data, f'{splitext(basename(system["psf"]))[0]}.dihedral_energy.html')

    # Write each conformation to a single PDB file with multiple frames
    # Load template PDB and its associated PSF
    pdb = PDB(system['pdb'])

    # Replace coordinates with the supplied ones and write out a new (temporary) PDB
    pdb_fname = f"{splitext(basename(system['psf']))[0]}.frames.pdb"
    pdb_fp = open(pdb_fname, 'w')

    for i in range(len(data)):
        coords = data[i]['coords']
        for j in range(len(pdb.atoms)):
            pdb.atoms[j].x = coords[j][0]
            pdb.atoms[j].y = coords[j][1]
            pdb.atoms[j].z = coords[j][2]
        pdb_fp.write('MODEL %d\n' % i)
        for energy_type in ('BOND', 'ANGLE', 'DIHED', 'IMPRP', 'ELECT', 'VDW'):
            pdb_fp.write(f'REMARK  42 {energy_type}: {data[i][energy_type]}\n')
        pdb.write(pdb_fp)
        pdb_fp.write('TER\nENDMDL\n')

    pdb_fp.close()
    print(f'Wrote coordinates of each frame to {pdb_fname}')

    # End goal here is to be able to dynamically edit the PRM file we are creating
    # for this ligand, and watch the MM energy surface change in relation to the QM
    # energy surface. Should in the end speed up this whole process.

if __name__ == '__main__':
    main()
