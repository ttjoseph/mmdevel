from AMBER import *
import sys
import argparse

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
        print >>sys.stderr, "Read %d atoms from %s." % (len(self.atoms), filename)
        
    def nuke_solvent(self):
        print >>sys.stderr, "Stripping solvent residues."
        self.atoms = [atom for atom in self.atoms if atom.resname not in self.solvent_residues]
    
    def write(self, fp):
        """Writes this PDB back out, in PDB format."""
        
        for atom in self.atoms:
            fp.write("ATOM  %5d %4s %3s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % \
                (atom.atomid, atom.atomname, atom.resname, atom.resid, \
                atom.x, atom.y, atom.z, \
                atom.occupancy, atom.tempfactor))
        


info("colorize_pdb: colorizes a PDB by residue")

parser = argparse.ArgumentParser()
parser.add_argument('pdb', help='PDB file to colorize')
parser.add_argument('value_file', help='At least one file of whitespace-separated numerical values - one output PDB will be generated per line', nargs='+', type=argparse.FileType('r'))
parser.add_argument('-n', '--normalization-factor', help='Normalization factor (all values are multiplied by this number)', type=float, default=1.0) 
args = parser.parse_args()

pdb = PDB(args.pdb)
pdb.nuke_solvent()
info("Multiplying all values by %.1f" % args.normalization_factor)

lines = []
for f in args.value_file:
    lines.extend(f.readlines())

for line in lines:
    try:
        values = [float(x) for x in line.split()]
    except ValueError:
        print "That be wack, yo."
        sys.exit()

    # Color the PDB by residue using occupancy and tempfactor fields
    # Fills in missing values with zero
    # We use a "virtual resid" because the resid might not start from 1, or be
    # contiguous, or make any sense at all.
    virtual_resid, last_resid = 0, pdb.atoms[0].resid
    for atom in pdb.atoms:
        if last_resid != atom.resid: virtual_resid += 1
        last_resid = atom.resid
        try:
            val = values[virtual_resid]
        except IndexError:
            val = 0
        val *= args.normalization_factor
        atom.occupancy = atom.tempfactor = round(val, 2)

    pdb.write(sys.stdout)
