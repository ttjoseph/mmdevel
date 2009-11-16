import os
import subprocess
import sys
import struct
import array
import gzip
import re
from time import sleep
from AMBER import *

def bomb(reason):
    '''Prints a message to stderr and then exits the program.'''
    print >>sys.stderr, "ERROR:", reason
    sys.exit(1)

def inclusive_range(i, j, stride=1):
    """Convenience wrapper for xrange to allow an inclusive range."""
    return xrange(i, j + 1, stride)

def zeros(length):
    '''Returns a list of zeros of the specified length.'''
    tmp = []
    for i in xrange(length):
        tmp.append(0)
    return tmp
    
def extract_interaction_energies_old(pairs, prmtop_filename, mdout_filenames):
    """Extracts interaction energies for specified residue pairs.
    Returns a dictionary of the form:

    {(a,b): [1.0, 2.0, -2.0, ...],
     (c,d): [1.0, 9.5, ...]}

    That is, with each tuple as a key and a list of energies as each
    associated value.
    """
    energies = {}
    for pair in pairs:
        energies[pair] = []
        
    prmtop = AmberSystem(prmtop_filename)
    
    # We make a new PairwiseEnergies object for each mdout file.
    # We may potentially have thousands of mdout files, especially if we
    # calculate energies every 0.2 ps like in that paper!
    # Supposedly at the end of the loop, Python will realize that
    # the object is no longer needed and garbage-collect it.
    for mdout in mdout_filenames:
        print "Working on %s..." % mdout
        en = PairwiseEnergies(prmtop.num_residues())
        # Shame we have to read the entire mdout file just to get a few values...
        en.add_from_mdout_file(mdout, 0, ['total'])
        for pair in pairs:
            # Grab the energy and save it
            energies[pair].append(en.get(Energies.TOTAL_RESIDUE_ENERGY, 'total', 0, pair[0], pair[1]))
    
    return energies       

def read_legacy_mmpbsa_parameters(filename):
    """Reads a MM_PBSA.pl parameter file and returns a dictionary with its parameters"""
    fp = faster_gzip_open(filename)
    category = ""
    params = {}
    for line in fp:
        # Strip comments and extraneous whitespace
        line = re.sub("#.+$", "", line).strip()
        # Length less than 3 means it can't possibly be a key-value pair
        if len(line) < 3:
            continue
        if line[0] == "@":
            category = line[1:3]
        else:
            (key, value) = re.split(" ", line)
            params[category + "_" + key] = value

    return params

def guess_num_atoms_from_prmtop(prmtop_filename):
    '''Reads the number of atoms in a system from its AMBER prmtop file.'''
    prmtop = faster_gzip_open(prmtop_filename)
    try:
        while True:
            line = prmtop.readline()
            if len(line) == 0: # EOF
                raise Exception("Couldn't guess the number of atoms from %s. Corrupt?" % prmtop_filename)
            if line[0:13] == "%FORMAT(10I8)":
                line = prmtop.readline()
                return int(line[0:8])
    finally:
        prmtop.close()

def guess_num_residues_from_prmtop(prmtop_filename):
    '''Reads the number of residues in a system from its AMBER prmtop file.'''
    prmtop = faster_gzip_open(prmtop_filename)
    try:
        while True:
            line = prmtop.readline()
            if len(line) == 0: # EOF
                raise Exception("Couldn't guess the number of residues from %s. Corrupt?" % prmtop_filename)
            if line[0:13] == "%FORMAT(10I8)":
                prmtop.readline()
                line = prmtop.readline()
                return int(line[8:16])
    finally:
        prmtop.close()

def guess_num_residues_from_mdout(mdout_filename):
    '''Reads the number of residues in a system from an output file (-o option 
    to sander) of an MD run.'''
    mdout = faster_gzip_open(mdout_filename)
    try:
        for line in mdout:
            if line[51:59] == "NRES   =":
                mdout.close()
                return int(line[59:67])
        raise Exception("Couldn't guess the number of residues from %s. Corrupt?" % mdout_filename)
    finally:
        mdout.close()

def guess_num_atoms_from_mdout(mdout_filename):
    '''Reads the number of atoms in a system from an output file (-o option 
    to sander) of an MD run.'''
    mdout = faster_gzip_open(mdout_filename)
    try:
        for line in mdout:
            if line[1:9] == "NATOM  =":
                return int(line[9:17])
        raise Exception("Couldn't guess the number of atoms from %s. Corrupt?" % mdout_filename)
    finally:
        mdout.close()

def navigate_to_prmtop_block(prmtop, block):
    search_str = "%%FLAG %s" % block
    while True:
        line = prmtop.readline()
        if len(line) == 0: # EOF
            raise Exception("Couldn't find the %s block in %s." % (block, prmtop_filename))
        # Is this the block we want?
        if line.find(search_str) == 0:
            return

def get_atom_names_from_prmtop(prmtop_filename):
    '''Returns a list of atom names from a prmtop file, as a string where each
    atom takes up 4 characters.
    
    A contrived example: "H   O   P   O1P O2PH5'1"'''

    num_atoms = guess_num_atoms_from_prmtop(prmtop_filename)
    return get_thing_names_from_prmtop(prmtop_filename, "ATOM_NAME", num_atoms)
    
def get_residue_names_from_prmtop(prmtop_filename):
    num_residues = guess_num_residues_from_prmtop(prmtop_filename)
    return get_thing_names_from_prmtop(prmtop_filename, "RESIDUE_LABEL", num_residues)

def get_thing_names_from_prmtop(prmtop_filename, thing, num_things):
    prmtop = faster_gzip_open(prmtop_filename)
    try:
        navigate_to_prmtop_block(prmtop, thing)
        prmtop.readline() # Eat format line
        thing_names = ""
        things_read = 0
        while things_read < num_things:
            line = prmtop.readline()
            # There are at most 20 things per line, so if we have fewer than
            # 20 things left to read, this must be a non-full line
            things_this_line = min(num_things - things_read, 20)
            thing_names += line[0 : 4 * things_this_line]
            things_read += things_this_line
        return thing_names
    finally:
        prmtop.close()
    raise Exception("Error finding %s names in %s" % (thing, prmtop_filename))

class Energies:
    TOTAL_RESIDUE_ENERGY = "TDC"
    BACKBONE_RESIDUE_ENERGY = "BDC"
    SIDECHAIN_RESIDUE_ENERGY = "SDC"

class PairwiseEnergies(Energies):
    '''Encapsulates residue-residue pairwise interaction energies.'''

    def add_from_mdout_file(self, mdout_filename, frame_id, \
        terms=['int', 'vdw', 'eel', 'pol', 'sa', 'total'], type=Energies.TOTAL_RESIDUE_ENERGY):
        mdout = faster_gzip_open(mdout_filename)
        # Python doesn't have native 2D arrays, and I don't feel like
        # adding a dependency to numpy
        energy_matrix = {}
        for term in terms:
            energy_matrix[term] = zeros(self.num_residues * self.num_residues)
        if 'total' not in terms:
            energy_matrix['total'] = zeros(self.num_residues * self.num_residues)

        try:
            energy = {}
            for line in mdout:
                # Total per residue energy includes backbone and sidechain contributions
                #if line[0:3] == "TDC" or line[0:3] == "BDC" or line[0:3] == "SDC":
                if line[0:3] == type:
                    resid1 = int(line[3:11])
                    resid2 = int(line[13:20])
                    energy['int'] = float(line[20:30])
                    energy['vdw'] = float(line[30:40])
                    energy['eel'] = float(line[40:50])
                    energy['pol'] = float(line[50:60])
                    energy['sa'] = float(line[60:70]) # This is actually an energy...
                    energy['total'] = 0
                    index = self.num_residues * (resid1 - 1) + (resid2 - 1)
                    # Save each term
                    for term in ['int', 'vdw', 'eel', 'pol', 'sa']:
                        if term in terms:
                            energy_matrix[term][index] = energy[term]
                        energy['total'] += energy[term]
                    # Save the total energy
                    energy_matrix['total'][index] = energy['total']
                
                    # I don't care about backbone vs sidechain and the various terms
                    # in the energy expression, so we've just summed them up.
                    # Otherwise this takes a silly amount of memory for a large complex.
        finally:
            mdout.close()

        # Actually save the data we read
        # term = 'total' # XXX: teh suck
        for term in terms:
            self.ensure_data_keys_exist(type, term, frame_id)
            self.data[type][term][frame_id] = energy_matrix[term]

    def __init__(self, num_residues):
        self.data = {}
        self.num_residues = num_residues
    
    def ensure_data_keys_exist(self, type, term, frame_id):
        """Ensure there is a place for this frame in self.data"""
        if type not in self.data:
            self.data[type] = {}
        if term not in self.data[type]:
            self.data[type][term] = {}
        if frame_id not in self.data[type][term]:
            # Make a 2D array, shiny and full of zeros
            self.data[type][term][frame_id] = zeros(self.num_residues * self.num_residues)
    
        
    def mean_over_all_frames(self, type, term, resid1, resid2):
        frame_list = self.data[type][term].keys()
        # term will probably be 'total'
        x = [self.data[type][term][frame_id][self.num_residues * (resid1 - 1) + (resid2 - 1)] \
            for frame_id in frame_list]
        return sum(x) / len(x)

    def set(self, value, type, term, frame_id, resid1, resid2):
        """Absurdly slow way to add an energy value."""
        self.ensure_data_keys_exist(type, term, frame_id)

        if resid1 > self.num_residues or resid2 > self.num_residues:
            print "Encountered a residue number which is greater than the" \
                + "specified maximum number of residues %d." % (self.num_residues)

        self.data[type][term][frame_id][self.num_residues * (resid1 - 1) + (resid2 - 1)] \
            = value
    
    def get(self, type, term, frame_id, resid1, resid2):
        return self.data[type][term][frame_id][self.num_residues * (resid1 - 1) + (resid2 - 1)]


class ResidueEnergies(Energies):
    '''Encapsulates a set of per-residue energies.'''
    def __init__(self):
        self.data = {}

    # energies.add(type='TDC', term='all', frame=1, resid=3, value=vdw_energy+...)
    def append(self, value, type, term, frame):
        '''Adds an energy value to the list of energy values.'''
        if type not in self.data:
            self.data[type] = {}
        if term not in self.data[type]:
            self.data[type][term] = {}
        if frame not in self.data[type][term]:
            self.data[type][term][frame] = []
        self.data[type][term][frame].append(value)

    def mean_over_all_frames(self, type, term, resid):
        count = 0
        total = 0
        for frame in self.data[type][term]:
            total += self.data[type][term][frame][resid]
            count += 1
        return total / count

class JobRunner:
    '''Primitive job queuing system for a single node with multiple processors.
    
    Use add() to add a SanderJob instance to the queue, and go() to start running
    everything.'''
    def __init__(self, max_cpu):
        self.max_cpu = max_cpu
        self.jobs = []

    def add(self, request):
        '''Adds a job to the queue.'''
        self.jobs.append(request)

    def go(self):
        '''Runs all the jobs (simultaneously as CPU limits allow)
        and waits until they complete. When this method returns, all the
        jobs are finished.'''
        queue = []

        # Prepare to run by adding jobs to the run queue
        for job in self.jobs:
            if job.prepare_for_launch() is False:
                raise Exception("There was an error preparing a job. Giving up.")
            queue.append(job)

        # Run the jobs
        running = []
        finished = []
        
        info("There are %d jobs to run, on up to %d processors. Fasten your seatbelt." % (len(queue), self.max_cpu))

        while len(queue) > 0:
            # Start as many jobs as we can
            while len(running) < self.max_cpu and len(queue) > 0:
                job = queue.pop(0)
                job.process = subprocess.Popen(job.args, stdout=subprocess.PIPE)
                # TODO: Check whether there was an error starting the job
                info("Started job: %s" % job.name)
                running.append(job)

            # Wait until one of them finishes
            zombie = []
            for job in running:
                job.process.poll()

                # If terminated, remove from run list and put in finished and
                # zombie lists (don't want to modify the run list while
                # iterating over it)
                if job.process.returncode is not None:
                    job.on_finished()
                    info("Done with job: %s" % job.name)
                    finished.append(job)
                    zombie.append(job)

            # Zombie list is just the jobs to remove from the run list
            for job in zombie:
                running.remove(job)

            # This is a busywait, so sleep a while
            sleep(1)


class SanderJob:
    '''Encapsulates a request to run a Sander job.'''
    
    @classmethod
    def molecular_mechanics(self, unique_id, prmtop_file, crd_file, footer=''):
        '''Factory method that creates a SanderJob meant for an MM calculation.'''
        x = SanderJob()
        x.args = [SanderJob.get_sander_path(), "-O"]
        x.name = "%s (MM)" % unique_id
        x.params = {'ntf': 1,
              'ntb': 0,
              'dielc': 1.0,
              'idecomp': 2,
              'igb': 2,
              'offset': 0.09,
              'intdiel': 1.0,
              'extdiel': 80.0,
              'gbsa': 2,
              'surften': 1.0,
              'cut': 999.0,
              'nsnb': 99999,
              'scnb': 2.0,
              'scee': 1.2,
              'imin': 1,
              'maxcyc': 1,
              'ncyc': 0}
        x.footer = footer
        x.frame_id = crd_file
        x.mdin_filename = os.path.join(os.path.dirname(crd_file), "%s.mdin" % unique_id)
        x.mdout_filename = os.path.join(os.path.dirname(crd_file), "%s.mdout" % unique_id)
        x.args.extend(["-i", x.mdin_filename, "-o", x.mdout_filename])
        x.args.extend(["-p", prmtop_file, "-c", crd_file])
        return x

    @classmethod
    def get_sander_path(cls):
        '''Returns the path of the sander executable, which should be $AMBERHOME/exe/sander
        
        Notably, we do not use the MPI version of sander. If we did, we would have to muck
        around with system-dependent ways to invoke mpirun or whatever.'''
        
        # Check AMBERHOME environment variable to find sander
        amberhome = os.environ.get("AMBERHOME")
        if amberhome is not None:
            executable = "%s/exe/sander" % amberhome
        
        if os.path.exists(executable) is False:
            raise Exception("Could not find sander executable. Did you set AMBERHOME?")
            
        return executable

    def prepare_for_launch(self):
        mdin = open(self.mdin_filename, "w")
        mdin.write(self.mdin_as_string())
        mdin.close()
        if os.path.exists(self.mdin_filename) is False:
            raise Exception("I tried to make a sander input file %s, but it didn't work." \
                % self.mdin_filename)
            return False
        return True

    def on_finished(self):
        ''' Called when the JobRunner finishes executing this job '''
        # gzip the output file
        os.system("gzip -f %s" % self.mdout_filename)
        self.mdout_filename = "%s.gz" % self.mdout_filename

    def set_param(self, k, v):
        '''Sets mdin parameters.'''
        self.params[k] = v
    
    def mdin_as_string(self):
        '''Returns a string containing the parameters in AMBER mdin format.'''
        s = '%s\n&cntrl\n' % self.name
        for k in sorted(self.params.keys()):
            s = s + " %s = %s,\n" % (k, self.params[k])
        s = s + "&end\n" + self.footer
        
        return s

class Trajectory:
    @classmethod
    def from_file(cls, filename):
        pass

# This class doesn't work yet...
class DCDTrajectory:
    INT_SIZE = 4
    
    def __init__(self, filename):
        self.filename = filename
        self.fp = faster_gzip_open(filename)
    
    def get(self, format):
        return struct.unpack(format, self.fp.read(struct.calcsize(format)))
    
    def read_header(self):
        magic, label = self.get("i4s")
        if label != "CORD":
            raise Exception("Something is fishy with DCD trajectory %s" % self.filename)
        icntrl = array.array('i')
        icntrl = self.get('21i')
        self.has_box = icntrl[10] == 1
        charmm_type = icntrl[11] == 1
        if self.has_box:
            info("There's a box in this DCD.")
        if charmm_type:
            info("This is a CHARMM-generated DCD.")
        self.num_frames = icntrl[0]
        start_timestep = icntrl[1]
        timestep_stride = icntrl[2]
        num_steps = icntrl[3]
        info("Frames: %s Start timestep: %d Stride: %d" % (self.num_frames, start_timestep, timestep_stride))
        junk, num_titles = self.get("2i")
        for i in xrange(0, num_titles):
            title, = self.get("80s")
            print title
        
        header2 = array.array('i')
        header2 = self.get('4i')
        self.num_atoms = header2[2]
        info("Number of atoms: %d" % self.num_atoms)
        # self.read_box()
        
    def read_box(self):
        if self.has_box:
            cell_array_size, = self.get('i') # should be 48 = 6 * sizeof(double)
            cell_array_size /= 8
            if cell_array_size != 6:
                raise Exception("Invalid cell array size %d" % cell_array_size)
            info("Cell array size: %d" % cell_array_size)
            cell = array.array('d')
            cell = self.get("%dd" % cell_array_size)
            info("Cell dimensions:", cell)
        
    def read_frame(self):
        self.read_box()
        x = array.array('f')
        y = array.array('f')
        z = array.array('f')
        junk, = self.get('i')
        x = self.get("%df" % self.num_atoms)
        junk, junk2 = self.get('2i')
        y = self.get("%df" % self.num_atoms)
        junk, junk2 = self.get('2i')
        z = self.get("%df" % self.num_atoms)
        junk, = self.get('i')
        return x, y, z

        