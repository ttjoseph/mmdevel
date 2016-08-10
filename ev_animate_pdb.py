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
        


info("ev_animate_pdb: animates PDB along eigenvector path")

if len(sys.argv) < 3:
    print >>sys.stderr, "Usage: <pdb> <value-filename>"
    sys.exit()


# Read in eigenvector - natom lines of 3-tuples
ev = []
for line in open(sys.argv[2]):
    try:
        values = [float(x) for x in line.strip().split()]
        ev.extend(values)
    except ValueError:
        print "That be wack, yo."
        sys.exit()

# Color the PDB by residue using occupancy and tempfactor fields
# Fills in missing values with zero
coeff = 0
step = 0.05 # 20 steps
scale = 150
while coeff <= 1:
    pdb = PDB(sys.argv[1])
    pdb.nuke_solvent()
    oldatoms = pdb.atoms
    i = 0
    for atom in pdb.atoms:
        pdb.atoms[i].x += ev[i*3] * coeff * scale
        pdb.atoms[i].y += ev[i*3+1] * coeff * scale
        pdb.atoms[i].x += ev[i*3+2] * coeff * scale
        i += 1

    pdb.write(sys.stdout)
    print "TER"
    print "END"
    coeff = coeff + step
