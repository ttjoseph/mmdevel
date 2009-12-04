from AMBER import *
import sys

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

if len(sys.argv) < 3:
    print >>sys.stderr, "Usage: <pdb> <value-filename>"
    sys.exit()

pdb = PDB(sys.argv[1])
pdb.nuke_solvent()


for line in open(sys.argv[2]):
    try:
        values = [float(x) for x in line.split()]
    except ValueError:
        print "That be wack, yo."
        sys.exit()

    # Color the PDB by residue using occupancy and tempfactor fields
    for atom in pdb.atoms:
        val = values[atom.resid - 1]
        # if val > 5: val = 5
        # if val < -1: val = -5
        val /= len(pdb.atoms); # BS normalization factor
        atom.occupancy = atom.tempfactor = round(val, 2)

    pdb.write(sys.stdout)