#!/usr/bin/env python
# Converts a PSF file in conjunction with CHARMM parameter (.prm) file into the binary
# prmtop format used by coulomb.c.
#
# Usage: <psf> <prm>
# Output is always in solute.top.tom
import array
import math
import sys
import struct
import re

BOND, ANGLE, DIHEDRAL = 1, 2, 4

def put_i32(f, data):
    """Writes a signed 32-bit integer to file handle f."""
    f.write(struct.pack('i', data))

def put_i32_array(f, data):
    # First, write number of elements
    put_i32(f, len(data))
    for x in data:
        put_i32(f, x)
    
def put_f32_array(f, data):
    # First, write number of elements
    put_i32(f, len(data))
    for x in data:
        f.write(struct.pack('f', x))

def put_u8_array(f, data):
    # First, write number of elements
    put_i32(f, len(data))
    for x in data:
        f.write(struct.pack('B', x))

def get_psf_block_header(f):
    """Formatted as I8 <string>"""
    l = f.readline()
    num = int(l[0:8].strip())
    kind = l[10:].strip()
    print "Encountered %s block with %d records" % (kind, num)
    return (num, kind)

def parse_ints(l):
    """Returns array of ints extracted from a single string. Each int should be
    8 characters, starting from first character in string."""
    
    num = len(l) - (len(l)%8)
    ret = []
    for i in xrange(0, num, 8):
        ret.append(int(l[i:i+8].strip()))
    return ret
    
def read_int_block(f, num_records, max_records_per_line):
    """Reads a block of integers from a PSF file.
    Assumes the file pointer is at the start of it."""
    data = []
    for index in xrange(num_records/max_records_per_line):
        l = f.readline()
        data.extend(parse_ints(l))

    if num_records % max_records_per_line != 0:
        l = f.readline()
        data.extend(parse_ints(l))
    
    return data
    
psf = open(sys.argv[1], "r")
prm = open(sys.argv[2], "r")
out = open("solute.top.tom", "w")

# Eat header
l = psf.readline()
if l[0:3] != "PSF":
    print >>sys.stderr, "%s doesn't look like a PSF file to me." % sys.argv[1]
    sys.exit(1)
psf.readline() # Blank line
(numlines, kind) = get_psf_block_header(psf)
for i in xrange(numlines): psf.readline()
psf.readline()

# !NATOM block contains, for each atom, its type and charge, as well as residue assignment
(num_atoms, kind) = get_psf_block_header(psf)
residue_map = []
charges = []
atom_types = []

# Parse each atom line
#          1         2         3         4         5         6
#0123456789012345678901234567890123456789012345678901234567890123456789
#       1 LIGH 1    ALA  N    NH3   -0.300000       14.0070           0
found_end_of_solute = False
num_solute_atoms = 0
last_resid = -1
num_residues = num_solute_residues = 0
cys_sg_list = [] # atom indices of CYS SG atoms
for index in xrange(num_atoms):
    l = psf.readline()
    atom_index = int(l[0:8].strip())
    segment = l[9:13].strip()
    resid = int(l[14:19].strip())
    resname = l[19:24].strip()
    atomname = l[24:29].strip()
    atomtype = l[29:34].strip()
    charge = float(l[35:44].strip())
    
    if resname == "CYS" and atomname == "SG":
        cys_sg_list.append(atom_index)
    
    # AMBER pre-multiplies by sqrt of coulomb constant
    # See http://ambermd.org/Questions/units.html
    charges.append(charge*math.sqrt(332.0636))
    atom_types.append(atomtype)
    
    # Are we in a new residue? If so, increment the residue count
    if last_resid != resid: num_residues += 1
    last_resid = resid
    residue_map.append(num_residues)
    
    if found_end_of_solute is False and resname == "TIP3":
        found_end_of_solute = True
        num_solute_residues = num_residues
        print "Looks like you have %d solute atoms (guessed by taking the atoms before the first water)." % num_solute_atoms
        
    if found_end_of_solute is False: num_solute_atoms += 1

psf.readline() # Blank line

# Read in the various bonds.
# !NBOND - 4 per line (8 numbers)
(num_bonds, kind) = get_psf_block_header(psf)
bonds = read_int_block(psf, num_bonds, 4)

psf.readline() # Blank line
# NTHETA block: angles, 3 per line, for 9 numbers
(num_angles, kind) = get_psf_block_header(psf)
angles = read_int_block(psf, num_angles, 3)

psf.readline() # Blank line
# NPHI block: dihedrals, 2 per line, for 8 numbers
(num_dihedrals, kind) = get_psf_block_header(psf)
dihedrals = read_int_block(psf, num_dihedrals, 2)

psf.readline() # Blank line
# IMPROPERS
(num_impropers, kind) = get_psf_block_header(psf)

if len(bonds)%2 != 0 or len(angles)%3 != 0:
    print >>sys.stderr, "I didn't parse bonds or angles properly. Whoops!"
    sys.exit(1)

# OK. Now we're done reading the PSF file.

# Generate list of distinct atom types
the_types = list(set(atom_types))
num_atom_types = len(the_types)
print "There are %d distinct atom types." % len(the_types)
# num_solute_atoms should be equal to the index of the first solvent atom

# We need to extract the LJ terms from the prm file. Happily, we don't
# need anything else from it (since we are only doing nonbonded energies).
l = ""
while l[0:4] != "NONB":
    l = prm.readline().strip()
if l[-1] == '-': l = prm.readline().strip() # Eat extra nonsense

nuke_comment_re = re.compile('!.*$') # Regexp to get rid of comments
# We assume the HBOND line is the end of the nonbonded parameters
prm_type, prm_epsilon, prm_rmin2 = [], [], []
prm_eps14, prm_rmin214 = [], []
while True:
    l = prm.readline()
    if l == "":
        break # Quit on EOF
    l = nuke_comment_re.sub('', l).strip()
    if l == "": continue
    if l[0:4] == "HBON": break
    atomtype = ""
    epsilon, rmin2, eps14, rmin214 = None, None, None, None
    try:
        (atomtype, dummy, epsilon, rmin2, dummy, eps14, rmin214) = l.split()
        eps14 = float(eps14.strip())
        rmin214 = float(rmin214.strip())
    except ValueError:
        (atomtype, dummy, epsilon, rmin2) = l.split()
    epsilon = float(epsilon.strip())
    rmin2 = float(rmin2.strip())
    prm_type.append(atomtype.strip())
    prm_epsilon.append(epsilon)
    prm_rmin2.append(rmin2)
    prm_eps14.append(eps14)
    prm_rmin214.append(rmin214)

# Now we have the LJ terms. AMBER does some extra preprocessing with these
# that CHARMM apparently doesn't. The LJ terms for each pair of atom types
# is precalculated before being stored in the prmtop file.
# Array of normal (non 1-4) LJ12 terms - one for each pair of atom types
# Array of normal LJ6 terms
nb_indices = [0] * num_atom_types**2
lj12 = [0] * ((num_atom_types**2 - num_atom_types)/2)
lj6 = [0] * ((num_atom_types**2 - num_atom_types)/2)
lj1214 = [0] * ((num_atom_types**2 - num_atom_types)/2)
lj614 = [0] * ((num_atom_types**2 - num_atom_types)/2)
counter = 0
for i in xrange(num_atom_types):
    for j in xrange(i):
        # prm_* variables contain nonbonded paramters for all atom types but
        # we will only store information for those types that occur in the structure
        idx_a = prm_type.index(the_types[i])
        idx_b = prm_type.index(the_types[j])
        eps_a, eps_b = prm_epsilon[idx_a], prm_epsilon[idx_b]
        rmin2_a, rmin2_b = prm_rmin2[idx_a], prm_rmin2[idx_b]
        eps14_a, eps14_b = prm_eps14[idx_a], prm_eps14[idx_b]
        rmin214_a, rmin214_b = prm_rmin214[idx_a], prm_rmin214[idx_b]
        # Figure out where in NBIndices this will be
        #     int nbparm_offs_i = Ntypes * (AtomTypeIndices[atom_i] - 1);
        #     int index = NBIndices[nbparm_offs_i+AtomTypeIndices[atom_j]-1] - 1;
        #     double thisEnergy = LJ12[index]*distRecip6*distRecip6 - LJ6[index]*distRecip6;
        # print "%s vs %s: %d %d" % (the_types[i], the_types[j], idx_a, idx_b)
        nb_indices[num_atom_types*i+j] = counter+1 # Don't forget stupid FORTRAN 1-based indexing
        nb_indices[num_atom_types*j+i] = counter+1
        # From CHARMM prm_all27:
        # V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
        # epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
        # Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
        eps = math.sqrt(eps_a * eps_b)
        rmin = rmin2_a + rmin2_b
        lj12[counter] = eps * rmin**12
        lj6[counter] = 2 * eps * rmin**6
        if None not in (eps14_a, eps14_b, rmin214_a, rmin214_b):
            eps14 = math.sqrt(eps14_a * eps14_b)
            rmin14 = rmin214_a + rmin214_b
            lj1214[counter] = eps14 * rmin14**12
            lj614[counter] = 2 * eps14 * rmin**6
        counter += 1
        
# Matrix of bond types
print "Making %d MB array for bond type cache. This is really slow." % (num_solute_atoms**2 / 1048576)
bond_type = array.array('B', [0] * num_solute_atoms * num_solute_atoms)
for i in xrange(0, len(bonds), 2):
    atom_i = bonds[i] - 1
    atom_j = bonds[i+1] - 1
    
    # if bonds[i] in cys_sg_list and bonds[i+1] in cys_sg_list:
    #     print "DISU atoms %d %d" % (bonds[i], bonds[i+1])
    
    if atom_i >= num_solute_atoms or atom_j >= num_solute_atoms:
        continue
    bond_type[num_solute_atoms*atom_i+atom_j] |= BOND
    bond_type[num_solute_atoms*atom_j+atom_i] |= BOND
    
for i in xrange(0, len(angles), 3):
    atom_i = angles[i] - 1
    atom_j = angles[i+2] - 1
    if atom_i >= num_solute_atoms or atom_j >= num_solute_atoms:
        continue
    bond_type[num_solute_atoms*atom_i+atom_j] |= ANGLE
    bond_type[num_solute_atoms*atom_j+atom_i] |= ANGLE

for i in xrange(0, len(dihedrals), 4):
    atom_i = dihedrals[i] - 1
    atom_j = dihedrals[i+3] - 1
    if atom_i >= num_solute_atoms or atom_j >= num_solute_atoms:
        continue
    # print "%d %d" % (atom_i, atom_j)
    bond_type[num_solute_atoms*atom_i+atom_j] |= DIHEDRAL
    bond_type[num_solute_atoms*atom_j+atom_i] |= DIHEDRAL

#
# Output format:
# Number of atoms, residues, atom types
put_i32(out, num_solute_atoms)
print "I suspect there are %d residues in the solute." % num_solute_residues
put_i32(out, num_solute_residues)
put_i32(out, num_atom_types)

# Output NBIndices
put_i32_array(out, nb_indices)

# Matrix of non-bonded terms indices - one cell for each atom type pair
#   This allows us to locate LJ params for a particular atom pair
# AtomTypeIndices: Array of atom type indices - atom types indexed by atom id
# AMBER does one-based indexing so that's what we had used before
atom_type_nums = [the_types.index(x)+1 for x in atom_types]
put_i32_array(out, atom_type_nums[0:num_solute_atoms])

# Output LJ12, LJ6
put_f32_array(out, lj12)
put_f32_array(out, lj6)

# Output Charges: Array of charges
put_f32_array(out, charges[0:num_solute_atoms])

# Output BondType
put_u8_array(out, bond_type)

# Output ResidueMap: Array of atom indexes for residues
put_i32_array(out, residue_map[0:num_solute_atoms])

# Extra Output: LJ12,1-4 and LJ6,1-4 parameters, which seems to be CHARMM's version
#   of the LJ 1-4 scaling (divide by 1.2) that AMBER uses
put_f32_array(out, lj1214)
put_f32_array(out, lj614)