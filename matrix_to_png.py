#!/usr/bin/env python
#
# Requires the pypng module, described at http://packages.python.org/pypng and
# available for download at https://github.com/drj11/pypng
#
# It just does a grayscale color ramp which is not that exciting.
import sys
import png

print >>sys.stderr, "matrix_to_png: Renders rescorrel matrix in text format to image\n"

if len(sys.argv) < 3:
    print >>sys.stderr, "Usage <rescorrel.txt> <rescorrel.png>"

f = open(sys.argv[1])
lines = f.readlines()
tmp = lines[0].split();
max_val = float(tmp[0])
min_val = max_val
for i in xrange(len(lines)):
    lines[i] = [float(n) for n in lines[i].split()]
    for n in lines[i]:
        if n > max_val: max_val = n
        if n < min_val: min_val = n
        
for i in xrange(len(lines)):
    for j in xrange(len(lines[0])):
        n = lines[i][j]
        n = n + min_val
        n = n / max_val
        # This makes zero values black and higher values more white
        lines[i][j] = int(n * 255)
    
print >>sys.stderr, "Output image dimensions:", len(lines), "by", len(lines[0]), "pixels"
    
image = png.from_array(lines, 'L')
image.save(sys.argv[2])
    
