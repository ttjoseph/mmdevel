#!/usr/bin/python
#
# Converts ptraj .pev essential dynamics to an easier to parse format, with
# one eigenvector per line.
import sys
import math

def load_frame(file, num_coords=None):
    # Eat comment line
    s = file.readline()
    if len(s) == 0:
        return None
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
    
    return frame

# Load average structure
avg = load_frame(sys.stdin)
num_coords = len(avg)

while True:
    frame = load_frame(sys.stdin, num_coords)
    if frame is None:
        break
    first = True
    for x in frame:
        if not first:
            sys.stdout.write(" ")
        sys.stdout.write(str(x))
        first = False
    sys.stdout.write("\n")