#!/usr/bin/python
#
# Converts ptraj .pev essential dynamics to an easier to parse format, with
# one eigenvector per line.
import sys
import math
import getopt

def load_frame(file, num_coords=None):
    # Eat comment line
    s = file.readline()
    if len(s) == 0:
        return (None, 0)
    dimension = file.readline().strip().split()
    dimension[0] = int(dimension[0])
    dimension[1] = float(dimension[1])
    
    # If number of coords wasn't specified, assume this is the average structure
    # frame and extract the number of coords
    if num_coords is None:
        num_coords = dimension[0]
    else:
        print >>sys.stderr, "Eigenvalue: %f" % dimension[1]
    
    # Up to seven coordinates per line
    num_lines = int(math.ceil(num_coords/7.0))
    frame = []
    # print "%d lines" % num_lines # DEBUG

    for i in xrange(num_lines):
        s = file.readline()
        # 6f11
        for j in xrange(0, 80, 11):
            x = s[j:j+11].strip()
            if len(x) == 0:
                break
            frame.append(float(x))
            # print "%f" % x
    
    return (frame, num_coords)

# N CA CB C
atom = ''
atom_order = ['N', 'CA', 'CB', 'C']
opts, args = getopt.getopt(sys.argv[1:], 'a:')
for o, a in opts:
    if o == '-a':
        if a not in atom_order:
            print >>sys.stderr, "Must choose one of N CA CB C, not %d" % a
            sys.exit(1)
        atom = a

# Load average structure - this guesses and sets num_coords
avg, num_coords = load_frame(sys.stdin)

# Delete CB for Gly residues (as prescribed by the missing_CB.dat file)
f = open("missing_CB.dat")
data = f.readlines()
missing_CB = [int(x) for x in data]

# Make a list of atom names corresponding to indices.
# Each group of 3 coords is xyz. So there should be numres*3*4 - numGly*3 coords.
num_residues = int((num_coords + 3*len(missing_CB))/4/3)
if atom != '':
    print >>sys.stderr, "Number of residues: %d" % num_residues
atom_list = []
for atom_id in xrange(1, num_residues*3+1):
    atom_list.append('N')
    atom_list.append('CA')
    if atom_id not in missing_CB:
        atom_list.append('CB')
    atom_list.append('C')

while True:
    frame, num_coords = load_frame(sys.stdin, num_coords)
    if frame is None:
        break
    first = True
    atom_id = 0
    for i in xrange(0, len(frame), 3):
        if not first:
            sys.stdout.write(" ")
        if atom == '' or atom_list[i/3] == atom:
            sys.stdout.write("%f %f %f" % (frame[i], frame[i+1], frame[i+2]))
        first = False
        atom_id += 1
    sys.stdout.write("\n")