#!/usr/bin/env python
from sys import stderr
print >>stderr, "energy_correl.py - Tom Joseph <thomas.joseph@mssm.edu>\n"
from AMBER import *
from MMPBSA import *
import math
import glob
import cPickle as pickle
import numpy

def print_percent_done(done_so_far, total, description = None):
    percent = int(float(done_so_far) / float(total) * 100)
    if done_so_far % (total / 20) == 0 or done_so_far == total:
        if description is not None:
            desc = " %s" % description
        else:
            desc = ""
        print >>sys.stderr, "%d%% done%s." % (percent, desc)

def load_energies_from_mdout(mdout_filename):
    fp = faster_gzip_open(mdout_filename)
    is_data = False
    row = 0
    for line in fp:
        # Sneakily extract the number of residues
        if is_data == False and line[51:59] == "NRES   =":
            num_residues = int(line[59:67])
            data = numpy.zeros((num_residues, num_residues))
        # Sentinel for start/end of pairwise interaction energies matrix
        if line[0:3] == "TOM":
            is_data = not is_data
        elif is_data:
            data[row,:] = [float(x) for x in line.split()]
            row += 1
    return data
            
def average_energies(mdout_filenames):
    """Returns a numpy 2D array with average values of pairwise residue interaction energies."""
    num_residues = guess_num_residues_from_mdout(mdout_filenames[0])
    total = numpy.zeros((num_residues, num_residues))
    for mdout_filename in mdout_filenames:
        total = total + load_energies_from_mdout(mdout_filename)
    return total / len(mdout_filenames)
    
def extract_interaction_energies(pairs, mdout_filenames):
    """Extracts energies only of those pairs for all mdout_filenames."""
    # Use len(pairs) by num_frames 2D array.
    # Each row is a pair, each column is a frame.
    energies = numpy.zeros((len(pairs), len(mdout_filenames)))
    mdout_counter = 0
    for mdout_filename in mdout_filenames:
        frame = load_energies_from_mdout(mdout_filename)
        pair_counter = 0
        for pair in pairs:
            energies[pair_counter,mdout_counter] = frame[pair[0], pair[1]]
            pair_counter += 1
        mdout_counter += 1
    return energies

def make_energies(mdout_filenames, cutoff=1.0):
    print >>sys.stderr, "Using %d mdout files." % len(mdout_filenames)
    # Find residue pairs satisfying a interaction energy threshold
    print >>sys.stderr, "Calculating average interaction energies."
    avg = average_energies(mdout_filenames)
    num_residues = numpy.shape(avg)[0]
    pairs = []
    for row in xrange(num_residues):
        for col in xrange(num_residues):
            # Ignore same and adjacent residue interactions
            if abs(row-col) > 1 and abs(avg[row,col]) > cutoff:
                # print "%d %d = %f" % (row, col, avg[row,col])
                pairs.append((row, col))
    print "Found %d pairs to examine, extracting energies." % len(pairs)
  
    # Extract selected interaction energies and do correlation analysis
    energies = extract_interaction_energies(pairs, mdout_filenames)
    print "Done extracting, calculating correlations."
    print "Shape of energies matrix:", energies.shape
    return (energies, pairs, num_residues)
    

def calculate_correlations(energies):
    """Returns a correlation matrix from a timeseries. Each pair has a row
    and each frame is a column."""
    num_pairs = energies.shape[0]
    num_frames = energies.shape[1]
    mean_energies = numpy.zeros(num_pairs)
    for pair in xrange(num_pairs):
        mean_energies[pair] = numpy.mean(energies[pair,:])
    
    correl = numpy.zeros((num_pairs, num_pairs))
    
    # Now, we are calculating correlations between pairs
    for ij in xrange(num_pairs):
        print_percent_done(ij + 1, num_pairs, "calculating correlations")
        for kl in xrange(num_pairs):
            if ij != kl:
                # Iterate over frames
                num, denom = 0, 0
                for t in xrange(num_frames):
                    a = energies[ij,t] - mean_energies[ij]
                    b = energies[kl,t] - mean_energies[kl]
                    num += a * b
                    denom += math.sqrt(a**2 * b**2)
                correl[ij,kl] = num / denom
    return correl

def trim_correl_matrix(pairs, correl, cutoff):
    """Removes all rows/columns from correlation matrix that don't have at least
    one element above the cutoff."""
    return (pairs, correl)
    
def project_correlations_onto_residues(pairs, num_residues, correl):
    """Projects correlations back onto residues, as in Kong and Karplus papers.
    
    RCij = sum_mn(C_mn x f(i,j))
    
    where f(i,j) == 1 iff residues i and j are involved in interaction (m,n))
    or (n,m).
    """
    residue_correl = numpy.zeros((num_residues, num_residues))
    
    # Iterate over residues
    for i in xrange(num_residues):
        print_percent_done(i, num_residues, "projecting correlation matrix")
        for j in xrange(num_residues):
            if i == j: continue
            # Iterate over correlations.
            # Rows and columns correspond to tuples in pairs variable
            for m in xrange(len(pairs)):
                for n in xrange(len(pairs)):
                    if m == n: continue
                    # Is i in interaction m and j in n, or i in n and j in m?
                    if (i in pairs[m] and j in pairs[n]) \
                        or (i in pairs[n] and j in pairs[m]):
                        residue_correl[i, j] += correl[m, n]
    return residue_correl
            
    
    
if __name__ == "__main__":
    mdout_filenames = glob.glob("snapshots/complex.crd.10??.mdout.gz")
    (energies, pairs, num_residues) = make_energies(mdout_filenames)
    correl = calculate_correlations(energies)
    print "Shape of correlation matrix:", correl.shape
    numpy.savetxt("correl.txt", correl, fmt="%11.6G ") # Extra space to ensure token separation
    # TODO: project correlation matrix back onto residues
    
