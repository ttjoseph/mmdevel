#!/usr/bin/env python
#
# Strategy:
# 
# Use SQLite to reorder the atoms so that atoms in the same residue are grouped together. We can then stuff this into
# AMBER's leap and make a prmtop out of it. Leap will further reorder the atoms, and we can make a new PDB.
# Now, given the old and new PDB files, we know the mapping of old to new, and can then reorder the trajectory (mdcrd format,
# converted from the DCD version using VMD).

import sys
import sqlite3 as sql
import gzip

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

        self.rename_opls_to_amber()

    def rename_opls_to_amber(self):
        '''Converts what I assume to be OPLS atom names to AMBER convention.
        
        This is done so that we can feed the resulting PDB into leap and not have it add heavy atoms or hydrogens,
        so we can generate a prmtop that corresponds to the original PDB. These atom names might not be anything
        special as they came from the Desmond toolset.'''
        if self.resname == 'ACE':
            if self.atomname == '1H': self.atomname = 'HH31'
            elif self.atomname == '2H': self.atomname = 'HH32'
            elif self.atomname == '3H': self.atomname = 'HH33'

        if self.resname == 'NMA':
            if self.atomname == '1HA': self.atomname = 'HH31'
            elif self.atomname == '2HA': self.atomname = 'HH32'
            elif self.atomname == '3HA': self.atomname = 'HH33'
            elif self.atomname == 'CA': self.atomname = 'CH3'
            self.resname = 'NME'

        if self.resname == 'GLU':
            if self.atomname == 'OXT': self.atomname = 'O'

        if self.resname == 'ASN':
            if self.atomname == 'H22': self.atomname = 'H'

        if self.resname == 'LYS':
            if self.atomname == 'H2': self.atomname = 'H'

        if self.resname == 'HIS':
            self.resname = 'HID'

class Molecule:
    '''Represents an entire molecule loaded from a PDB file.'''
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

    def keep_only_atomname(self, atomname):
        '''Keeps only a certain type of atom. (e.g. CA, CB, N, ...)'''
        print >>sys.stderr, "Keeping only %s atoms." % atomname
        self.atoms = [atom for atom in self.atoms if atom.atomname == atomname]
        
def atom_record(atomid, atomname, resname, chain, resid, x, y, z, occupancy, tempfactor):
    """Returns a PDB ATOM record string ('card' in FORTRAN parlance)"""
    return "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f                  " % \
        (atomid, atomname, resname, chain, resid, x, y, z, occupancy, tempfactor)
    
def print_molecule_from_sql(conn):
    '''Given an sqlite connection, dumps the molecule in table 'mol' to stdout in PDB format.'''
    c = conn.cursor()
    c.execute("SELECT * FROM mol ORDER BY chain, resid, atomid")
    conn.commit()
    recs = c.fetchall()
    for i in xrange(len(recs)):
        rec = recs[i]
        (atomid, atomname, resname, chain, resid, x, y, z, occupancy, tempfactor) = rec
        print atom_record(i+1, atomname, resname, chain, resid, x, y, z, occupancy, tempfactor)

def load_mdcrd_frame(mdcrd, num_atoms, box=True):
    """Loads a frame from an open mdcrd file, returning a list of coordinates that is 3*num_atoms in length.
    
    Box information is discarded."""
    
    num_lines = num_atoms * 3 / 10
    if num_atoms*3%10 != 0: num_lines += 1
    
    frame = []
    for i in xrange(num_lines):
        s = mdcrd.readline()
        # If we've (possibly prematurely) reached EOF, this is not a complete frame, so return None
        if s == '': return None
        frame += [float(x) for x in s.split()]
    # Discard box information if it is present
    if box: mdcrd.readline()
    
    if len(frame) != num_atoms * 3:
        print >>sys.stderr, "There are %d atoms but %d coordinates read from the mdcrd. Aw." % (num_atoms, len(frame))
    
    return frame
    
def print_frame(frame):
    '''Prints a trajectory snapshot in mdcrd format.'''
    for i in xrange(len(frame)):
        if i != 0 and i % 10 == 0: print ""
        print "%8.3f" % frame[i],
        

if __name__ == '__main__':
    ordered_mol, scrambled_mdcrd = None, None
    scrambled_mol = Molecule(sys.argv[1])
    scrambled_mol.nuke_solvent()
    
    # If we just want to reorder a molecule, there will be only one command-line argument.
    # If there are second and third arguments, we want to map the atoms to each other and convert a trajectory.
    if len(sys.argv) == 4:
        ordered_mol = Molecule(sys.argv[2])
        ordered_mol.nuke_solvent()
        if len(scrambled_mol.atoms) != len(ordered_mol.atoms):
            print >>sys.stderr, "Number of atoms, %d vs %d does not match, fool." % (len(scrambled_mol.atoms), len(ordered_mol.atoms))
            sys.exit(1)
        # Open the trajectory file and discard the header line
        scrambled_mdcrd = gzip.open(sys.argv[3], 'r')
        scrambled_mdcrd.readline()
    
    # Create a temporary in-memory table to store the structure
    conn = sql.connect(":memory:")
    c = conn.cursor()
    c.execute("CREATE TABLE mol (atomid integer, atomname text, resname text, chain text, resid integer, x real, y real, z real, occupancy real, tempfactor real)")

    # Insert atoms into in-memory table
    insert_query = "INSERT INTO mol (atomid, atomname, resname, chain, resid, x, y, z, occupancy, tempfactor) VALUES (:atomid, :atomname, :resname, :chain, :resid, :x, :y, :z, :occupancy, :tempfactor)"
    for a in scrambled_mol.atoms:
        c.execute(insert_query, {"atomid": a.atomid, "atomname": a.atomname, "resname": a.resname, "chain": a.chain, "resid": a.resid, "x": a.x, "y": a.y, "z": a.z, "occupancy": a.occupancy, "tempfactor": a.tempfactor})
    conn.commit()
    
    # If we only wanted to reorder the PDB, just do that and exit
    if ordered_mol is None:
        print_molecule_from_sql(conn)
        sys.exit(0)
        
    # We must now create a 1:1 mapping between scrambled_mol and ordered_mol.
    # The trajectory is scrambled, so for each atom of a frame, we want to do:
    #   Ordered[i] = Scrambled[Map[i]]
    # where i is the atom id in the ordered trajectory.
    
    # XXX:
    # Running a SELECT for each atom in ordered_mol is incredibly stupid. I admit this.
    # A much faster and cleaner solution would be to do a JOIN. But I am lazy.
    scrambled_to_ordered = {}
    for i in xrange(len(ordered_mol.atoms)):
        # Try to find this atom in the scrambled atom table
        a = ordered_mol.atoms[i]
        q = "SELECT * FROM mol WHERE atomname = :atomname AND resname = :resname AND chain = :chain AND resid = :resid"
        c.execute(q, {"atomname": a.atomname, "resname": a.resname, "chain": a.chain, "resid": a.resid})
        conn.commit()
        recs = c.fetchall()
        if len(recs) > 1:
            print >>sys.stderr, "More than one atom found in scrambled_mol for a given ordered atom!"
            sys.exit(1)
        r = recs[0]
        # print "atomid: ordered[%d] = scrambled[%d]" % (a.atomid, r[0])
        scrambled_to_ordered[r[0]] = a.atomid
        
    print >>sys.stderr, "Done creating mapping between scrambled and ordered molecules."
    
    # Now we can iterate through our scrambled trajectory and our mapping tells us what maps to what
    num_atoms = len(scrambled_mol.atoms)
    print "Generated by %s" % sys.argv[0]
    while True:
        scrambled_frame = load_mdcrd_frame(scrambled_mdcrd, num_atoms)
        ordered_frame = [0] * num_atoms * 3
        if scrambled_frame is None: break
        for i in xrange(num_atoms):
            # atom ids are 1-based. I hate fortran.
            ordered_i = scrambled_to_ordered[i+1]-1
            ordered_frame[ordered_i*3] = scrambled_frame[i*3]
            ordered_frame[ordered_i*3+1] = scrambled_frame[i*3+1]
            ordered_frame[ordered_i*3+2] = scrambled_frame[i*3+2]
        print_frame(ordered_frame)
