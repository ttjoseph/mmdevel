#!/usr/bin/env python
#
# Python library for molecular modeling stuff
import sys
from math import pi
import numpy as np

BOHRS_PER_ANGSTROM = 1.88972612
KCAL_MOL_PER_HARTREE = 627.509474

class AtomRecord:
    '''Represents a single atom.'''    
    def __init__(self, line):
        # We are treating atomid as some arbitrary string, because it might be > 99999,
        # and some software deals with this by putting hex or possibly other stuff in
        # this field.
        self.atomid = line[6:11].strip()
        self.atomname = line[12:16].strip()
        # Allow 4-character residue names, which are quite common
        self.resname = line[17:21].strip()
        self.chain = line[21]
        #if self.chain == " ":
        #     self.chain = "X"
        # We are specifically leaving resid as a string because in the real world it 
        # is not guaranteed to be a parseable into an integer (e.g. VMD atomselect writepdb)
        self.resid = line[22:26].strip()
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = float(line[54:60])
        self.tempfactor = float(line[60:66])
        self.segid = line[72:76].strip()
        
    def write(self, f):
        """Writes this atom to an ATOM record to an open file."""
        line = f'ATOM  {self.atomid:>5s} {self.atomname:<4s} {self.resname:4s}{self.chain:1s}{self.resid:>4s}    {self.x:8.3f}{self.y:8.3f}{self.z:8.3f}{self.occupancy:6.3f}{self.tempfactor:6.3f}      {self.segid:4s}\n'
        # write() expects bytes, not a string, so we need to specify the encoding
        f.write(line.encode('utf-8'))

class ShadyPDB:
    """Shady representation of a Protein Data Bank (PDB) file that might not conform to the standard.

    This was initially made for regenerate_psf.

    Why do we have our own code for parsing PDB files? Aren't there a ton of libraries?

    The PDB file format is outdated and generally crap, because it uses fixed-width fields,
    as a legacy of punch-card style I/O from Fortran. VMD copes with this by using hex in the
    residue ID field. This means that we can no longer assume there are decimal integers only
    in this field. Both MDAnalysis and Bio.PDB, to their credit, appear to do a pretty good
    job of following the PDB spec, but this means we can't use them when there are more than
    9999 residues in a single PDB segment.
    """
    def __init__(self, filename):
        self.atoms = []
        self.load_from_pdb(filename)

    def create_indexes(self):
        """Build a mapping of residue index (NOT resid, which is some arbitrary string)
        to atom index (NOT atomid, which is also some arbitrary string), and segids
        to residue indexes.

        Ultimately we are doing everything by index instead of resid or atomid because
        we can't trust that resids are either integers or unique outside a single
        segment. Same for atomid.
        """
        self.segid_to_resindex = {}
        self.resindex_to_atomindex = []
        last_resid, last_segid, last_resindex = self.atoms[0].resid, self.atoms[0].segid, 0
        this_residue_atomidx = []
        
        for a_idx in range(len(self.atoms)):
            a_resid = self.atoms[a_idx].resid
            a_segid = self.atoms[a_idx].segid
            # Is this residue different from the last one?
            # Can happen when resid changes or when segid changes (resid might be the same,
            # for example when there is a one-residue segment)
            if a_resid != last_resid or a_segid != last_segid:
                # The residue (and possibly segment) changed, so save the current indexing we've done
                if last_segid not in self.segid_to_resindex:
                    self.segid_to_resindex[last_segid] = []
                self.segid_to_resindex[last_segid].append(last_resindex)
                self.resindex_to_atomindex.append(this_residue_atomidx)

                last_resindex += 1
                last_resid = a_resid
                last_segid = a_segid
                this_residue_atomidx = []

            this_residue_atomidx.append(a_idx)
        # We won't trigger these actions on the last atom (or even once, if there's only one atom) so we do them now
        self.segid_to_resindex[last_segid].append(last_resindex)
        self.resindex_to_atomindex.append(this_residue_atomidx)

    def set_segid(self, atomindexes, new_segid):
        """Sets the segid for a specified subset of atoms."""
        for a_idx in atomindexes:
            self.atoms[a_idx].segid = new_segid

    def renumber_subset(self, atomindexes, renumber_residues=True, renumber_atoms=True):
        """Renumber a subset of atoms, starting from 1."""
        new_atomid, new_resid = 1, 1
        sorted_atomindexes = list(sorted(atomindexes))
        last_resid = self.atoms[sorted_atomindexes[0]].resid
        for a_idx in sorted_atomindexes:
            if self.atoms[a_idx].resid != last_resid:
                last_resid = self.atoms[a_idx].resid
                new_resid += 1
            self.atoms[a_idx].atomid = str(new_atomid)
            self.atoms[a_idx].resid = str(new_resid)
            new_atomid += 1

    def write_to_pdb(self, file, atomindexes=None, flush=True):
        we_opened_file = False
        if type(file) == 'str':
            f = open(file, 'w')
            we_opened_file = True
        else:
            f = file

        if atomindexes is None:
            atomindexes = range(len(self.atoms))
        
        for a_idx in atomindexes:
            self.atoms[a_idx].write(f)

        # Ensure everything is written before we return.
        # Could be that the caller is invoking VMD or something, and we want to ensure all is written.
        if flush is True:
            f.flush()

        if we_opened_file:
            f.close()
        
    def load_from_pdb(self, filename):
        coords_set = set()
        pdb = open(filename)
        lines = pdb.readlines()
        for line in lines:
            if line.startswith("ATOM  "):
                atom = AtomRecord(line)
                self.atoms.append(atom)
                coords_set.add((atom.x, atom.y, atom.z))
        pdb.close()
        self.create_indexes()
        print(f"ShadyPDB: Read {len(self.atoms)} atoms from {filename}.", file=sys.stderr)
        if len(coords_set) != len(self.atoms):
            print(f"ShadyPDB: Warning: More than one atom in {filename} has the same coordinates! This will ruin your simulation!")


class CoorVel(object):
    """NAMD-format binary .coor or .vel file.

    These have the same format: number of atoms (int32), followed by xyz coordinates for each atom as doubles.
    Importantly, we assume everything is little endian, because most machines running NAMD these days are
    in fact little endian.
    """

    def __init__(self):
        self.num_atoms = np.uint32(0)
        self.coords = np.zeros(0)

    def blank(self, num_atoms):
        self.num_atoms = np.uint32(num_atoms)
        self.coords = np.zeros(num_atoms*3, dtype=np.dtype('<f8'))

    def load(self, fname):
        """Loads .coor or .vel file into this object."""
        # So format is totatoms (32-bit integer) then 3 doubles for each atom
        # Assume little endian. Not sure if NAMD etc always write little endian or what.
        num_atoms = np.fromfile(fname, dtype=np.dtype('<i4'), count=1, offset=0)
        self.num_atoms = num_atoms[0]
        self.coords = np.fromfile(fname, dtype=np.dtype('<f8'), offset=4, count=self.num_atoms*3)

    def write(self, fname):
        """Writes this object in binary format suitable for NAMD."""
        f = open(fname, 'wb')
        f.write(self.num_atoms.astype(np.dtype('<i4')).tobytes())
        f.write(self.coords.tobytes())
        f.close()

    def excise(self, atom_indices):
        """Removes specified atoms in place."""
        # Convert atom indices to coordinate indices
        # Ex: [2, 3] -> [2*3, 2*3+1, 2*3+2, 3*3, 3*3+1, 3*3+2]
        x = atom_indices * 3
        y = atom_indices * 3 + 1
        z = atom_indices * 3 + 2
        xyz = np.concatenate((x, y, z), axis=None)
        xyz.sort()
        # print(f'Deleting coord indices: {xyz}', file=sys.stderr)
        # print(self.coords[xyz].reshape((-1, 3)), file=sys.stderr)
        self.coords = np.delete(self.coords, xyz)
        self.num_atoms -= len(atom_indices)


def pbc_charge_correction(ion_charge, ion_radius, box_length, dielectric=80, verbose=False):
    """Correction, in kcal/mol, for the charging of an ion in periodic boundary conditions.

    From Eqn 17 in Simonson and Roux, Molecular Simulation, 2016
    https://www.tandfonline.com/doi/full/10.1080/08927022.2015.1121544

    Of note, the first term of the correction includes a division by the dielectrc constant.
    Not all authors have done this, particularly for simulations with nonpolarizable force fields.
    My naive suspicion is that polarizable force fields provide dielectric screening that
    nonpolarizable force fields do not, and therefore the calculated ion charging energies
    are lower. So the correction needs to take this into account. Since no (new) dipoles are
    induced with a nonpolarizable force field, each periodic image point charge experiences the
    full interaction energy with the others.
    """
    ZETA = -2.837297 # units are Hartree * Bohr
    # 1 Angstrom is 1.88972612 Bohr, so we will convert sizes to Bohr
    box_length *= BOHRS_PER_ANGSTROM
    ion_radius *= BOHRS_PER_ANGSTROM
    term1 = (ion_charge**2 * np.abs(ZETA))/(2*dielectric*box_length)
    if ion_radius > 0:
        term2 = (2*pi * ion_charge**2 * ion_radius**2) / (3 * box_length**3)
        term2 *= (dielectric - 1) / dielectric
    else:
        if verbose:
            print('pbc_charge_correction: No ion radius supplied, so not doing second term', file=sys.stderr)
        term2 = 0
    # term3 = O(R**2 L**-3 dielectric**-3)
    # 1 Bohr is 0.529177211 Angstrom
    # 1 Hartree, for one particle -> 627.509474 kcal/mol
    # This correction is in Hartrees, but we want kcal/mol
    if verbose:
        print(f'term1 = {term1*KCAL_MOL_PER_HARTREE:.3f} kcal/mol; term2 = {term2*KCAL_MOL_PER_HARTREE:.3f} kcal/mol',
            file=sys.stderr)
    correction = term1 + term2
    correction *= KCAL_MOL_PER_HARTREE
    return correction

