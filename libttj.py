#!/usr/bin/env python
#
# Python library for molecular modeling stuff
import sys
import os
from glob import glob
from collections import defaultdict
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

        # This field is not part of the PDB standard.
        # We need a residue index value because the PDB resid is not guaranteed to be
        # unique! Some software just dumps in 
        self.resindex = None
        
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
        self.segid_to_resindex = defaultdict(list)
        self.resindex_to_atomindex = []
        last_resid, last_segid, last_resindex = self.atoms[0].resid, self.atoms[0].segid, 0
        this_residue_atomidx, this_residue_atomnames = [], []
        
        for a_idx in range(len(self.atoms)):
            a_resid = self.atoms[a_idx].resid
            a_segid = self.atoms[a_idx].segid
            # Keep track of the atom names we've seen so far in this residue
            a_atomname = self.atoms[a_idx].atomname
            # Is this residue different from the last one?
            # Can happen when resid changes or when segid changes (resid might be the same,
            # for example when there is a one-residue segment)
            # or when we see the same atom name again
            if a_resid != last_resid or a_segid != last_segid or a_atomname in this_residue_atomnames:
                # The residue (and possibly segment) changed, so save the current indexing we've done
                if last_segid not in self.segid_to_resindex:
                    self.segid_to_resindex[last_segid] = []
                self.segid_to_resindex[last_segid].append(last_resindex)
                self.resindex_to_atomindex.append(this_residue_atomidx)
                last_resindex += 1
                last_resid = a_resid
                last_segid = a_segid
                this_residue_atomidx, this_residue_atomnames = [], []

            # Save this atom's resindex (which is confusingly named and separate from resid).
            # We do it here because we might have rolled over resindex in the if block just above.
            self.atoms[a_idx].resindex = last_resindex

            # Also save this atom's index within the atoms array, so the caller has a single
            # unambiguous way to identify it
            self.atoms[a_idx].atomindex = a_idx

            this_residue_atomidx.append(a_idx)
            this_residue_atomnames.append(a_atomname)
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
        last_resindex = self.atoms[sorted_atomindexes[0]].resindex
        
        for a_idx in sorted_atomindexes:
            a = self.atoms[a_idx]
            if self.atoms[a_idx].resindex != last_resindex:
                last_resindex = self.atoms[a_idx].resindex
                new_resid += 1
            if renumber_atoms:
                self.atoms[a_idx].atomid = str(new_atomid)
            if renumber_residues:
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
                if (atom.x, atom.y, atom.z) in coords_set:
                    print(f"ShadyPDB: Coordinates of {atom.atomname} {atom.atomid} are duplicated. This is really bad!", file=sys.stderr)
                coords_set.add((atom.x, atom.y, atom.z))
        pdb.close()
        self.create_indexes()
        print(f"ShadyPDB: Read {len(self.atoms)} atoms from {filename}.", file=sys.stderr)
        if len(coords_set) != len(self.atoms):
            print(f"ShadyPDB: Warning: More than one atom in {filename} has the same coordinates! This will ruin your simulation!",
                file=sys.stderr)

    def set_coords_from_coorvel(self, coorvel):
        """Overwrites our xyz coordinates with those from a CoorVel object that was hopefully loaded from a .coor file.

        This is useful for remap_coorvel.py, where we may want to use the atom names from the starting PDB but the
        coordinates from the in-progress simulation to decide what is the atom correspondence.
        """
        if len(coorvel.coords) != (len(self.atoms) * 3):
            print(f'ShadyPDB: The PDB has {len(self.atoms)} atoms but the coorvel object has coords for {len(coorvel.coords)/3} atoms!',
                file=sys.stderr)
            return
        
        for i in range(len(self.atoms)):
            self.atoms[i].x, self.atoms[i].y, self.atoms[i].z = coorvel.coords[i*3], coorvel.coords[i*3+1], coorvel.coords[i*3+2]


class ShadyPSF(object):
    """Shady representation of a NAMD PSF file.
    
    Not necessarily CHARMM or X-PLOR variants, though those might work by accident.
    """
    # Atom tuple indexes, so we don't need to resort to magic numbers
    ATOMID = 0
    SEGID = 1
    RESID = 2
    RESNAME = 3
    ATOMNAME = 4
    ATOMTYPE = 5
    CHARGE = 6
    WEIGHT = 7

    def __init__(self, filename):
        self.atoms = []
        self.load_from_psf(filename)

    def load_from_psf(self, filename):
        psf = open(filename, 'r')

        # Eat header
        l = psf.readline()
        if l[0:3] != "PSF":
            print("%s doesn't look like a PSF file to me." % sys.argv[1], file=sys.stderr)
            sys.exit(1)
        
        # Parse each block
        while True:
            records, kind = self.parse_block(psf)
            if kind is None:
                break
            elif kind.startswith('!NTITL'):
                pass
            elif kind.startswith('!NATOM'):
                self.atoms = []
                for record in records:
                    atomid = record[ShadyPSF.ATOMID]
                    segid = record[ShadyPSF.SEGID]
                    resid = record[ShadyPSF.RESID]
                    resname = record[ShadyPSF.RESNAME]
                    atomname = record[ShadyPSF.ATOMNAME]
                    atomtype = record[ShadyPSF.ATOMTYPE]
                    charge = float(record[ShadyPSF.CHARGE])
                    weight = float(record[ShadyPSF.WEIGHT])
                    self.atoms.append(record)
            elif kind.startswith('!NBOND'):
                self.bonds = ShadyPSF.separate_into_tuples(records, 2)
            elif kind.startswith('!NTHET'):
                self.angles = ShadyPSF.separate_into_tuples(records, 3)
            elif kind.startswith('!NPHI'):
                self.dihedrals = ShadyPSF.separate_into_tuples(records, 4)
            elif kind.startswith('!NIMPHI'):
                self.impropers = ShadyPSF.separate_into_tuples(records, 4)
            elif kind.startswith('!NCRTERM'):
                self.crossterms = ShadyPSF.separate_into_tuples(records, 8)
            elif kind.startswith('!NDON'):
                self.donors = ShadyPSF.separate_into_tuples(records, 2)
            elif kind.startswith('!NACC'):
                self.acceptors = ShadyPSF.separate_into_tuples(records, 2)
            elif kind.startswith('!NNB'):
                # It's not really one thing per tuple, but this will allow us to
                # spit it back out properly
                self.exclusions = ShadyPSF.separate_into_tuples(records, 1)
            else:
                print(f'ShadyPSF: Unknown PSF block kind {kind}, so I am discarding it. You should look into this.',
                    file=sys.stderr)

            # print(f'ShadyPSF: {filename}: {kind} block has {len(records)} records.', file=sys.stderr)

    @classmethod
    def separate_into_tuples(cls, records, tuple_size):
        raw_data = []
        for record in records:
            raw_data.extend([int(s) for s in record])
        if len(raw_data) % tuple_size != 0:
            print(f'ShadyPSF: Error: Number of entries {len(raw_data)} is not a multiple of {tuple_size}, but it should be.',
                file=sys.stderr)
            sys.exit(1)
        
        data = []
        for i in range(0, len(raw_data), 3):
            data.append(tuple(raw_data[i:i+tuple_size]))
        return data

    def parse_block(self, f):
        """Parse a PSF block from an open file handle.
        
        Returns a list of lists of tokens. It's up to the caller to parse further than that.
        """
        num_things, kind = ShadyPSF.get_psf_block_header(f)
        if kind is None:
            return [], None

        records = []
        # We don't know how many lines there are a priori
        # If an ATOM or TITLE block, it's one record per line
        # If any of the other kinds of block, it's a variable number of records per line
        if kind.startswith(('!NATOM', '!NTITL')):
            for _ in range(num_things):
                line = f.readline()
                records.append(line.strip().split())
        else:
            num_things_left = num_things
            # The number of things is not the number of tokens
            # e.g. One angle is actually three tokens
            if kind.startswith('!NBOND'):
                num_things_left *= 2
            elif kind.startswith('!NTHET'):
                num_things_left *= 3
            elif kind.startswith(('!NPHI', '!NIMPHI')):
                num_things_left *= 4
            elif kind.startswith('!NCRTERM'):
                num_things_left *= 8
            elif kind.startswith('!NNB'):
                # Exclusions list
                # Number of tokens is num_things_left + number of atoms
                if len(self.atoms) == 0:
                    print(f'ShadyPSF: Trying to parse NNB (exclusions) block but have not parsed atoms yet.',
                        file=sys.stderr)
                    print(f'ShadyPSF: This means the PSF file is probably corrupt.', file=sys.stderr)
                    sys.exit(1)
                num_things_left += len(self.atoms)
            while num_things_left > 0:
                line = f.readline()
                tokens = line.strip().split()
                num_things_left -= len(tokens)
                records.append(tokens)

        return records, kind

    @classmethod
    def get_psf_block_header(cls, f):
        """Parse PSF block header which is formatted as I8 <string>.
        
        This is not an instance method because it doesn't need to be.
        """
        # Eat blank lines or nonsense data until we get to a block header, and give up on EOF
        while True:
            l = f.readline()
            # Give up at EOF
            if l == '':
                break
            tokens = l.strip().split()
            if len(tokens) >= 2 and tokens[1][0] == '!':
                num, kind = int(tokens[0]), tokens[1].strip(':')
                print("ShadyPSF: Encountered %s block with %d records" % (kind, num), file=sys.stderr)
                return num, kind
        return 0, None

    @classmethod
    def read_int_block(cls, f, num_records, max_records_per_line):
        """Reads a block of integers from a PSF file.
        Assumes the file pointer is at the start of it.
        
        This is not an instance method because it doesn't need to be.
        """
        data = []
        for index in range(num_records/max_records_per_line):
            l = f.readline()
            data.extend(ShadyPSF.parse_ints(l))

        if num_records % max_records_per_line != 0:
            l = f.readline()
            data.extend(ShadyPSF.parse_ints(l))
        
        return data

    @classmethod
    def parse_ints(cls, l):
        """Returns array of ints extracted from a single string. Each int should be
        8 characters, starting from first character in string."""
        
        num = len(l) - (len(l)%8)
        ret = []
        for i in range(0, num, 8):
            ret.append(int(l[i:i+8].strip()))
        return ret



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


def parse_fepout(fnames, verbose=False):
    """Given a list of .fepout filenames (generated by NAMD), return lists of cumulative and delta ∆G as well as
    the lambda values used.
    """
    fepenergy, deltas, lambdas = list(), list(), list()
    lines = list()

    if isinstance(fnames, str):
        fnames = [fnames,]

    for fname in fnames:
        if verbose:
            print(f"parse_fepout: Processing {fname}...", file=sys.stderr)
        if os.path.exists(fname) is False:
            print(f"parse_fepout: I was asked to parse file {fname} but it doesn't seem to exist",
                file=sys.stderr)
            return None, None, None
        if os.path.getsize(fname) == 0:
            print(f"parse_fepout: File {fname} exists but has zero size", file=sys.stderr)
            return None, None, None        
        with open(fname) as f:
            theselines = f.readlines()
            goodlines = []
            line_counter = 1
            for line in theselines:
                # Check for malformed FepEnergy: lines
                if line.startswith('FepEnergy'):
                    tokens = line.split()
                    if len(tokens) < 10:
                        print(f"parse_fepout: Messed up FepEnergy line in {fname}:{line_counter} with line length {len(line)}", file=sys.stderr)
                        print(f'Here are the first 200 characters in the line ({len(tokens)} tokens):', file=sys.stderr)
                        print(line[:200], file=sys.stderr)
                        continue
                goodlines.append(line)

            lines.extend(goodlines)

    for line in lines:
        tokens = line.split()
        if line.startswith('FepEnergy'):
            # A single sample
            # The energy is always in the tenth column
            if len(tokens) < 10:
                print(f"parse_fepout: Skipping messed up FepEnergy line", file=sys.stderr)
                continue
            energy = float(tokens[9])
            if energy > 20:
                print(f"parse_fepout: Skipping messed up FepEnergy energy {tokens[9]}", file=sys.stderr)
                continue
            fepenergy.append(energy)
        elif line.startswith('#Free'):
            # delta-G at end of window, in the twelfth column
            deltas.append(float(tokens[11]))
        elif line.startswith('#NEW FEP WINDOW'):
            # Choose the smallest value of the various reported lambdas to represent this window
            lambda1, lambda2 = float(tokens[6]), float(tokens[8])
            lambdaIDWS = float(tokens[10]) if 'IDWS' in line else np.Inf
            # lambdas.append(np.min([lambda1, lambda2, lambdaIDWS]))
            lambdas.append(lambda1)

    if len(lambdas) != len(deltas):
        print(f"parse_fepout: While processing {fnames}:")
        print(f"parse_fepout: We have {len(lambdas)} lambdas but {len(deltas)} delta-G values", file=sys.stderr)
        return None, None, None

    # Sort each list by the lambdas, since there's no guarantee the lambdas are in any
    # particular order
    fepenergy = np.array(fepenergy)
    deltas = np.array(deltas)
    lambdas = np.array(lambdas)
    sorted_indices = np.argsort(lambdas)

    return fepenergy[sorted_indices], deltas[sorted_indices], lambdas[sorted_indices]


# Return the first matching glob of files
# This will allow us to use (for example) either dis5A000_fwd.fepout, or dis5A???.fepout
def try_globs(*globspecs):
    for globspec in globspecs:
        fnames = glob(globspec)
        if len(fnames) > 0:
            break
    return sorted(fnames)


def find_fepouts(dirnames, fepdirname, prefixes, verbose=False):
    good_dirnames = set()
    fwd_fnames, bwd_fnames = defaultdict(list), defaultdict(list)
    for dirname in dirnames:
        if os.path.exists(dirname) is False or os.path.isdir(dirname) is False:
            print('Could not find directory {}.'.format(dirname), file=sys.stderr)
            continue

        # Try each prefix (dis5A, dis5B, ...)
        for prefix in prefixes.split(','):
            fwd_fname = try_globs(f'{dirname}/{fepdirname}/{prefix}???_fwd.fepout',
                f'{dirname}/{fepdirname}/{prefix}???.fepout')
            bwd_fname = try_globs(f'{dirname}/{fepdirname}/{prefix}???_bwd.fepout')

            # For each directory, there are two dicts: one for forward fepouts, and another for backward
            # fepouts. We require at least a forward fepout. If no backward fepouts, just save an empty list.
            # This ensures that both dicts have the same number of lists of fepouts, so we can iterate
            # through them later properly and not lose which backward fepout sets correspond to which
            # forward fepout sets.
            if len(fwd_fname) == 0:
                # print('Could not find forward fepouts in {}'.format(dirname))
                continue
            else:
                if verbose:
                    print(f"Found forward fepouts with prefix {prefix} in {dirname}", file=sys.stderr)
                fwd_fnames[dirname].append(fwd_fname)
            if len(bwd_fname) == 0:
                if verbose:
                    print('Could not find backward fepouts in {}'.format(dirname), file=sys.stderr)
                bwd_fnames[dirname].append([])
            else:
                bwd_fnames[dirname].append(bwd_fname)

            # Only care about this directory if we found a fepout
            good_dirnames.add(dirname)

    return list(sorted(good_dirnames)), fwd_fnames, bwd_fnames