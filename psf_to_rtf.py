#!/usr/bin/env python
#
# Generates a CHARMM residue topology from a PSF file
import MDAnalysis as mda
import argparse

def main():
    ap = argparse.ArgumentParser(description='Generate an RTF from a PSF')
    ap.add_argument('psf')
    args = ap.parse_args()

    u = mda.Universe(args.psf)
    # Print atoms and bonds. This is straightforward enough.
    print("""* Converted from:
* %s
*
36 1

RESI %s      0.000
GROUP
""" % (args.psf, u.atoms[0].resname))
    for a in u.atoms:
        print('ATOM %-6s %-6s %8.6f' % (a.name, a.type, a.charge))
    print('')
    for b in u.bonds:
        print('BOND %-4s %-4s' % (b.atoms[0].name, b.atoms[1].name))

    # Not so sure if MDAnalysis gives you PSF CMAP terms, so let's just read the things
    # directly ourselves
    cmap_tokens = []
    total_num_cmap_tokens = 0
    with open(args.psf) as f:
        in_cmap_block = False
        
        for line in f:
            tokens = line.strip().split()
            if len(tokens) > 1 and '!NCRTERM' in line:
                in_cmap_block = True
                total_num_cmap_tokens = int(tokens[0]) * 8
            elif len(tokens) > 1 and tokens[1].startswith('!'):
                in_cmap_block = False
            elif in_cmap_block and len(cmap_tokens) < total_num_cmap_tokens:
                cmap_tokens.extend([int(x) for x in tokens])
    print('')
    for i in range(0, len(cmap_tokens), 8):
        s = 'CMAP'
        for j in range(i+0, i+8):
            # PSFs start counting atoms at 1
            s += ' ' + u.atoms[cmap_tokens[j] - 1].name
        print(s)
        print('')
        print('END')

if __name__ == '__main__':
    import sys
    sys.exit(int(main() or 0))