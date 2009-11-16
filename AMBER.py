import gzip
import re
import sys
import os
import subprocess

def info(s):
    print >>sys.stderr, "INFO: %s" % s

def faster_gzip_open(filename):
    """Python 2.4's gzip module is absurdly slow: an order of magnitude slower
    than the command-line gzip. So we use the latter instead. Supposedly newer Pythons
    have a faster gzip, so we use it if we are using Python 2.5 or higher.

    Also, automatically checks for filename.gz and opens that if it exists, instead.
    """

    is_gzipped = False

    if filename[-3:] == ".gz":
        is_gzipped = True
    elif os.path.exists("%s.gz" % filename):
        info("Opening %s.gz instead of %s." % (filename, filename))
        filename = "%s.gz" % filename
        is_gzipped = True

    if is_gzipped:
        if sys.version_info[0] == 2:
            if sys.version_info[1] < 5:
                return subprocess.Popen("gzip -dc %s" % filename, shell=True, stdout=subprocess.PIPE).stdout
            else:
                return gzip.open(filename, "r")
        else:
            raise Exception("faster_gzip_open: I am too stupid to work with your Python.")            
    else:
        return open(filename, "r")

class AmberSystem:
    # POINTERS block indices, stolen from the prmtop format spec on ambermd.org
    NATOM  = 0 # total number of atoms 
    NTYPES = 1 # total number of distinct atom types
    NBONH  = 2 # number of bonds containing hydrogen
    MBONA  = 3 # number of bonds not containing hydrogen
    NTHETH = 4 # number of angles containing hydrogen
    MTHETA = 5 # number of angles not containing hydrogen
    NPHIH  = 6 # number of dihedrals containing hydrogen
    MPHIA  = 7 # number of dihedrals not containing hydrogen
    NHPARM = 8 # currently not used
    NPARM  = 9 # currently not used
    NEXT   = 10 # number of excluded atoms
    NRES   = 11 # number of residues
    NBONA  = 12 # MBONA + number of constraint bonds
    NTHETA = 13 # MTHETA + number of constraint angles
    NPHIA  = 14 # MPHIA + number of constraint dihedrals
    NUMBND = 15 # number of unique bond types
    NUMANG = 16 # number of unique angle types
    NPTRA  = 17 # number of unique dihedral types
    NATYP  = 18 # number of atom types in parameter file, see SOLTY below
    NPHB   = 19 # number of distinct 10-12 hydrogen bond pair types
    IFPERT = 20 # set to 1 if perturbation info is to be read in
    NBPER  = 21 # number of bonds to be perturbed
    NGPER  = 22 # number of angles to be perturbed
    NDPER  = 23 # number of dihedrals to be perturbed
    MBPER  = 24 # number of bonds with atoms completely in perturbed group
    MGPER  = 25 # number of angles with atoms completely in perturbed group
    MDPER  = 26 # number of dihedrals with atoms completely in perturbed groups
    IFBOX  = 27 # set to 1 if standard periodic box, 2 when truncated octahedral
    NMXRS  = 28 # number of atoms in the largest residue
    IFCAP  = 29 # set to 1 if the CAP option from edit was specified
    
    def __init__(self, prmtop_filename, rst_filename=None):
        """Initializes this instance with AMBER prmtop and rst files."""
        self.blocks = {}
        self.block_list = []
        self.formats = {}
        self.format_strings = {}
        self.x, self.y, self.z = [], [], []
        self.box = None
        self.name = prmtop_filename
        
        self.load_prmtop(prmtop_filename)
        if rst_filename is not None:
            self.load_rst(rst_filename)
                
    def load_prmtop(self, filename):
        """Loads an AMBER prmtop file into this instance."""
        
        # These regular expressions are used in parsing FLAG and FORMAT lines
        block_name_re = re.compile('%FLAG\s+(\w+)')
        self.format_re = re.compile('%FORMAT\((\d+)(.)([\d\.]+)')
        
        # Read header line, TODO check version of file
        fp = faster_gzip_open(filename)
        self.header = fp.readline()
        
        # Read first line first
        line = fp.readline()
        while True:
            if line == "": break
            name = block_name_re.match(line).groups()[0]
            line = fp.readline()
            assert(line != "", "%s truncated after %%FLAG card." % filename)
            format = self.format_re.match(line).groups()
            tokens_per_line = int(format[0])
            token_type = format[1]
            token_length = float(format[2])
            self.formats[name] = (tokens_per_line, token_type, token_length)
            self.format_strings[name] = line
            self.new_block(name, line)
            
            # Read lines until we get one that starts with %, which signifies
            # a new block
            line = fp.readline()
            if line == "": break
            while line[0] is not "%":
                # Subtract 1 from length to ignore newline
                for i in xrange(0, len(line) - 1, int(token_length)):
                    self.blocks[name].append(line[i:i + int(token_length)])
                line = fp.readline()
                if line == "": break
        
            # Converts tokens to the correct type and saves them in arrays
            # I: integer; a: alphanumeric; E: float
            if token_type == "I":
                self.blocks[name] = [int(x) for x in self.blocks[name]]
            elif token_type == "E":
                self.blocks[name] = [float(x) for x in self.blocks[name]]
        
    def load_rst(self, filename):
        """Loads an AMBER restart file into this instance."""
        
        assert 'POINTERS' in self.blocks, "POINTERS block missing - valid prmtop not loaded"
        
        fp = faster_gzip_open(filename)
        fp.readline() # Eat header line
        
        coords_left = self.blocks['POINTERS'][AmberSystem.NATOM] * 3
        assert(coords_left == int(fp.readline()) * 3, \
            "Number of atoms %d in %s differs from that specified in prmtop" \
                % (coords_left, filename))
        
        # There are 6 coordinates, 12 characters each per line
        while coords_left > 0:
            line = fp.readline()
            self.x.append(float(line[0:12]))
            self.y.append(float(line[12:24]))
            self.z.append(float(line[24:36]))
            coords_left -= 3
            if coords_left >= 3:
                self.x.append(float(line[36:48]))
                self.y.append(float(line[48:60]))
                self.z.append(float(line[60:72]))
                coords_left -= 3
                
        assert(coords_left == 0, \
            "Number of coordinates read is not a multiple of 3: is %s corrupt?" % \
            filename)
        
        # If a box is specified in the prmtop, read it in
        if self.blocks['POINTERS'][AmberSystem.IFBOX] != 0:
            self.box = []
            line = fp.readline()
            self.box.append(float(line[0:12]))
            self.box.append(float(line[12:24]))
            self.box.append(float(line[24:36]))
            self.box.append(float(line[36:48]))
            self.box.append(float(line[48:60]))
            self.box.append(float(line[60:72]))
        
    def save_prmtop(self, filename):
        """Saves the prmtop part of this AmberSystem into a prmtop file."""
        
        fp = open(filename, "w")
        fp.write(self.header)
        
        # Eases converting from Fortran to Python/C format specifiers
        type_table = {'I': 'd', 'E': 'E', 'a': 's'}

        # The following ugly code deals with printing tokens in the correct format
        for name in self.block_list:
            fp.write("%%FLAG %s\n" % name)
            fp.write(self.format_strings[name])
            format = self.format_re.match(self.format_strings[name]).groups()
            tokens_per_line = int(format[0])
            token_type = type_table[format[1]]
            token_length = float(format[2])
            token_length_int = int(token_length)
            # Print x.y vs x length specifiers correctly
            if float(token_length_int) == token_length:
                token_length_dec = ""
            else:
                token_length_dec = "%.1f" % (token_length - token_length_int)
                token_length_dec = token_length_dec[1:]
            s = "%%%d%s%s" % (token_length_int, token_length_dec, token_type)

            count = 0
            for x in self.blocks[name]:
                fp.write(s % x)
                count += 1
                if count % tokens_per_line == 0:
                    fp.write("\n")
            
            if count % tokens_per_line != 0:
                fp.write("\n")
        
        fp.close()
        
    def save_rst(self, filename):
        """Saves the coordinates and box information (if present) to a restart file."""
        
        assert len(self.x) == len(self.y) == len(self.z)

        fp = open(filename, "w")
        # Write header
        fp.write("Processed by softcore_setup.py\n")
        fp.write("%6d\n" % len(self.x))

        # Write coordinates
        for i in xrange(len(self.x)):
            fp.write("%12.7f" % self.x[i])
            fp.write("%12.7f" % self.y[i])
            fp.write("%12.7f" % self.z[i])
            if i % 2 == 1: fp.write("\n")
                
        if len(self.x) % 2 == 1: fp.write("\n")
            
        if self.box is not None:
            for n in self.box:
                fp.write("%12.7f" % n)
            fp.write("\n")
        
        fp.close()
    
    def save_pdb(self, sclist=None, filename=None):
        """Makes a PDB that LEaP will read back in without complaint. This
        differs from ambpdb/ambmask in that it writes beta and occupancy fields
        that are set to 1.0 if that atom is in sclist, and 0.0 if not.
        
        This means:
        - Don't mangle names (e.g. don't do 1H5' instead of H5'1)
        - Use single quotes instead of asterisks
          (which can be accomplished simply by not mangling the prmtop names)
        - Insert a TER record after every MOLECULE
        - Have an END record at the end, I guess
        """
        
        if filename is not None:
            fp = open(filename, "w")
        else:
            fp = sys.stdout
        
        # Iterate over atoms and determine which residue
        # it is in. Save this information in arrays.
        residues = []
        residue_counter = self.num_residues() - 1
        for atom in xrange(self.num_atoms(), 0, -1):
            residues.append(residue_counter + 1)
            if atom == self.blocks['RESIDUE_POINTER'][residue_counter]:
                residue_counter -= 1
        residues.reverse()
        
        # We count down the number of atoms in each molecule as we go so we
        # can print TER records as needed
        if 'ATOMS_PER_MOLECULE' in self.blocks:
            have_molecules = True
            current_molecule = 0
            atoms_left = self.blocks['ATOMS_PER_MOLECULE'][current_molecule]
        else:
            have_molecules = False
        
        # Print the ATOM records
        for atom in xrange(self.num_atoms()):
            atom_name = self.blocks['ATOM_NAME'][atom]
            res_name = self.blocks['RESIDUE_LABEL'][residues[atom] - 1]
            
            if sclist is not None and (atom + 1) in sclist:
                is_softcore = 1.0
            else:
                is_softcore = 0.0
            
            fp.write("ATOM  %5d %4s %3s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % \
                (atom + 1, atom_name, res_name, residues[atom], \
                self.x[atom], self.y[atom], self.z[atom], \
                is_softcore, is_softcore))
            
            # If there's a ATOMS_PER_MOLECULE block we know where to print the TERs    
            if have_molecules:
                atoms_left -= 1
                if atoms_left == 0:
                    fp.write("TER\n")
                    current_molecule += 1
                    if current_molecule < len(self.blocks['ATOMS_PER_MOLECULE']):
                        atoms_left = self.blocks['ATOMS_PER_MOLECULE'][current_molecule]
        fp.write("END\n")
        
        if filename is not None:
            fp.close()
        
    # The hash we are using to store the blocks does not maintain
    # the order of the keys. We keep an ordered list
    # of keys so when we save the prmtop the blocks will be in that order.
    def new_block(self, name, format_string):
        """Creates a new prmtop block, setting up ancillary data structures as needed."""
        
        if name not in self.blocks:
            assert name not in self.block_list, \
                "AmberSystem block_list inconsistent with blocks"
            self.block_list.append(name)
            self.blocks[name] = []
            self.format_strings[name] = format_string
            
    def residue_for_atom(self, atom):
        """Returns the residue number of a given atom index, or None if not 
        found."""
        
        assert('RESIDUE_POINTER' in self.blocks, \
            "No RESIDUE_POINTER block - Corrupt prmtop")
        if atom < 1 or atom > self.num_atoms(): return None
        
        for residue in xrange(len(self.blocks['RESIDUE_POINTER']), 0, -1):
            if atom >= self.blocks['RESIDUE_POINTER'][residue - 1]:
                return residue
        
        raise Error("Bizarre residue boundaries - corrupt prmtop?")
        
    def residues_for_atoms(self, atoms):
        """Returns a hash of residue numbers to the provided atom numbers."""
        
        tmp = {}
        for atom in atoms:
            residue = self.residue_for_atom(atom)
            if residue not in tmp:
                tmp[residue] = [atom]
            else:
                tmp[residue].append(atom)
        return tmp
        
    def make_ambmask(self, atoms):
        """Converts a list of atoms into a more human-friendly ambmask string."""
        
        atoms_by_residue = self.residues_for_atoms(atoms)
        mask = ""
        first_residue = True
        for residue in atoms_by_residue:
            if not first_residue: mask += " | "
            first_residue = False
            mask += ":%d@" % residue
            atom_names = [self.blocks['ATOM_NAME'][i - 1].strip() for i in \
                atoms_by_residue[residue]]
            first_atom = True
            for atom_name in sorted(atom_names):
                if not first_atom: mask += ","
                first_atom = False
                mask += atom_name
        return mask
    
    def num_atoms(self):
        """Returns the number of atoms. Must have loaded a prmtop (with
        load_prmtop)."""

        assert('POINTERS' in self.blocks, "Corrupt prmtop")
        return int(self.blocks['POINTERS'][AmberSystem.NATOM])
        
    def num_residues(self):
        """Returns the number of residues."""
        
        assert('POINTERS' in self.blocks, "Corrupt prmtop")
        return int(self.blocks['POINTERS'][AmberSystem.NRES])
        
    def load_trajectory(self, filename):
        self.trajectory = AmberTrajectory.from_crd_file(filename, self.num_atoms(), \
            box is not None)
        

class AmberTrajectory:
    '''Encapsulates an AMBER trajectory.'''

    @classmethod
    def from_frames(cls, frames):
        ''' Constructs an AMBER trajectory object from a list of Frames. '''
        t = AmberTrajectory()
        t.atoms_per_frame = frames[0].num_atoms()
        t.has_box = frames[0].has_box
        t.frames = frames
        return t

    @classmethod
    def from_crd_file(cls, filename, atoms_per_frame, has_box = True):
        ''' Reads an AMBER trajectory from a file. '''
        t = AmberTrajectory()
        t.filename = filename
        t.filetype = "crd"
        t.atoms_per_frame = atoms_per_frame
        t.has_box = has_box
        return t

    def __init__(self):
        # Have to initialize this stuff here and not directly under the
        # class declaration because in that case they become class
        # variables, not instance variables, and hilarity ensues
        self.wanted_frames = []
        self.frames = []

    def get_frames(self, wanted_frames = [], wanted_atoms = []):
        '''Gets specified frames, containing specified atoms.'''
        frames = self.load_frames_from_file(wanted_frames)

        if wanted_atoms == []:
            return frames

        return AmberTrajectory.from_frames([x.extract_atoms(wanted_atoms) for x in frames])
        
    def snapshotify(self, prefix):
        frame_id = 1
        fp = faster_gzip_open(self.filename)
        # Eat header line
        fp.readline()
        
        while True:
            frame = Frame(fp = fp, atoms_per_frame = self.atoms_per_frame)
            if frame.num_atoms() == 0: # Happens when we hit EOF
                break
            
            info("Writing frame %d, with length %d" % (frame_id, len(frame.x)))
            filename = "%s.crd.%d" % (prefix, frame_id)
            out = open(filename, "w")
            # Header line
            out.write("Sippin on gin and juice\n")
            out.write("%d\n" % frame.num_atoms())
            # f12.7, six coords per line
            max_even_number = frame.num_atoms() - frame.num_atoms() % 2
            for i in xrange(0, max_even_number, 2):
                out.write("%12.7f%12.7f%12.7f" % (frame.x[i], frame.y[i], frame.z[i]))
                out.write("%12.7f%12.7f%12.7f\n" % (frame.x[i + 1], frame.y[i + 1], frame.z[i + 1]))

            # If the number of atoms is odd, we missed the last one above,
            # so write that one too
            if frame.num_atoms() % 2 != 0:
                i = frame.num_atoms() - 1
                out.write("%12.7f%12.7f%12.7f" % (frame.x[i], frame.y[i], frame.z[i]))

            frame_id += 1

    def load_frames_from_file(self, wanted_frames=[]):
        ''' Reads selected frames from a trajectory file. '''
        # Return cached frames if we can
        if self.frames != [] and self.wanted_frames == wanted_frames:
            info("Hey, I have those frames already! No need to read them again!")
            return self.frames

        # TODO: actually look at file to detect whether its gzipped?
        fp = faster_gzip_open(self.filename)

        # Eat header line
        fp.readline()
        # Read frames, keeping only the ones we want
        frame_count = 0
        while True:
            # Frame constructor will read in a frame from the file handle
            frame = Frame(fp = fp, atoms_per_frame = self.atoms_per_frame, has_box=self.has_box)
            if frame.num_atoms() == 0: # Happens when we hit EOF
                break
            if len(wanted_frames) == 0 or frame_count in wanted_frames:
                frame.id = frame_count
		#import pdb; pdb.set_trace()
                self.frames.append(frame)
                #info("Keeping frame %d." % (frame_count))
            frame_count = frame_count + 1
            if frame_count > max(wanted_frames):
                break

        fp.close()
        self.wanted_frames = wanted_frames
        return self.frames

    def write_crd_files(self, prefix):
        '''Writes the frames in this trajectory to separate crd files
        {prefix}.crd.{i}'''
        frame_id = 1 # People are used to counting from 1 in this situation
        for frame in self.frames:
            info("Writing frame %d, with length %d" % (frame_id, len(frame.x)))
            filename = "%s.crd.%d" % (prefix, frame_id)
            out = open(filename, "w")
            # Header line
            out.write("Welcome to flavor country\n")
            out.write("%d\n" % frame.num_atoms())
            # f12.7, six coords per line
            max_even_number = frame.num_atoms() - frame.num_atoms() % 2
            for i in xrange(0, max_even_number, 2):
                out.write("%12.7f%12.7f%12.7f" % (frame.x[i], frame.y[i], frame.z[i]))
                out.write("%12.7f%12.7f%12.7f\n" % (frame.x[i + 1], frame.y[i + 1], frame.z[i + 1]))

            # If the number of atoms is odd, we missed the last one above,
            # so write that one too
            if frame.num_atoms() % 2 != 0:
                i = frame.num_atoms() - 1
                out.write("%12.7f%12.7f%12.7f" % (frame.x[i], frame.y[i], frame.z[i]))

            frame_id += 1


class Frame:
    '''Encapsulates a trajectory frame.'''

    def init(self):
        self.x = []
        self.y = []
        self.z = []
        self.box_a = 0
        self.box_b = 0
        self.box_c = 0
        self.atoms_per_frame = 0
        self.id = 0

    def __init__(self, atoms_per_frame, fp=None, has_box=True):
        self.init()
        self.atoms_per_frame = atoms_per_frame
        self.has_box = has_box
        if fp != None:
            self.from_file(fp)

    def extract_atoms(self, indices=[]):
        '''Returns a new Frame containing only the specifed atoms.'''

        # If no atom indices are specified, assume that caller wants them all
        if len(indices) == 0:
            indices = range(0, self.num_atoms())

        # Indices will start from 1, but we start from 0
        # my_indices = sorted([x - 1 for x in indices])

        # Make a new Frame that just contains the requested atoms
        #print "extract_atoms: Atoms: %d, indices: %d" % (len(self.x), len(indices))
        
        frame = Frame(atoms_per_frame = self.atoms_per_frame)
        frame.x = [self.x[i] for i in indices]
        frame.y = [self.y[i] for i in indices]
        frame.z = [self.z[i] for i in indices]
        return frame

    def from_file(self, fp):
        '''Reads a frame from the trajectory file pointed to by filehandle fp.'''
        num_tokens = 0
        coords = []
        tokens_per_frame = self.atoms_per_frame * 3

        while True:
            line = fp.readline()

            if len(line) == 0: # Must have reached end of file
                x = y = z = [] # This frame is not valid, so get rid of coords
                return

            for column in range(0, 80, 8):
                try:
                    coords.append(float(line[column : column + 7]))
                    num_tokens = num_tokens + 1
                except ValueError:
                    # The line is not full of tokens...
                    pass

            if num_tokens >= tokens_per_frame:
                if self.has_box:
                    line = fp.readline()
                    (self.box_a, self.box_b, self.box_c) = line.split()
                break

        self.x = [coords[i] for i in range(len(coords)) if i % 3 == 0]
        self.y = [coords[i] for i in range(len(coords)) if i % 3 == 1]
        self.z = [coords[i] for i in range(len(coords)) if i % 3 == 2]


    def num_atoms(self):
        ''' Returns the number of atoms in this frame. '''
        return len(self.x)

    def sanity_check(self):
        if len(self.x) != len(self.y) or len(self.y) != len(self.z):
            return False
        return True


def ambmask_to_atom_list(prmtop, rst, ambmask):
    """Runs $AMBERHOME/bin/ambmask, extracts atom indices from its output, and
    returns them as a list"""

    executable = "%s/bin/ambmask" % os.environ['AMBERHOME']
    if os.path.exists(executable) is False:
        raise Error("ambmask executable not found. I checked" \
            " %s. Did you set the AMBERHOME environment variable correctly?" \
            % executable)
    p = subprocess.Popen([executable, '-p', prmtop, '-c', rst, '-find', ambmask],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = p.communicate()[0]
    if "Error" in out:
        raise Error("There's an error in ambmask %s:\n\n%s" % (ambmask, out))
    # Return the atom indices only
    out = out.split('\n')
    return [int(line[6:11]) for line in out if line[0:6] == 'ATOM  ']
