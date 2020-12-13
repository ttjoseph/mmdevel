#!/usr/bin/env python
#
# With the help of VMD/psfgen, takes a PDB that is probably from an existing simulation but maybe slightly
# modified and makes a new PSF to go with it. Importantly, the PDB must have already been in suitable shape
# for a simulation.
#
# You might ask: Why not use AutoPSF? Because AutoPSF is unbelievably slow for this use case.
#
# Tom Joseph, University of Pennsylvania

import sys
from os import remove
import argparse
import pathlib
import tempfile
import subprocess
import MDAnalysis as mda


def main():
    ap = argparse.ArgumentParser(description='Remakes a PSF/PDB from existing PDB. Helpful if you want to avoid the CHARMM-GUI rigamarole')
    ap.add_argument('pdb', help='PDB file to use')
    ap.add_argument('toppar_txt', help='File containing newline-separated list of psfgen-acceptable topology files')
    ap.add_argument('out_prefix', help='Prefix for output PSF and PDB files')
    ap.add_argument('--chunk-size', type=int, default=9999, help='Maximum number of residues in a segment that you think psfgen or whatever can handle')
    ap.add_argument('--show-vmd-output', action='store_true', help='Show the output of VMD, which may be long')
    ap.add_argument('--show-psfgen-script', action='store_true', help='Show the psfgen script that will be passed to VMD')
    ap.add_argument('--vmdbin', default='vmd', help='VMD binary')
    ap.add_argument('--keep-temp-pdbs', action='store_true', help='Keep the temporary PDBs made for psfgen so you can look at them')
    args = ap.parse_args()

    # psfgen needs to have input PDBs chopped up into segments and regurgitated into it, much like a helpless baby bird.
    # I can't figure out how to get VMD to deal with a number of residues > 9999 when we start with a PDB so we will do the segment
    # splitting here, in Python. (This often is necessary when there is a bunch of water in a simulation system.)
    out_seg_suffixes = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    u = mda.Universe(args.pdb)
    segs_to_process = []
    for seg_i in range(len(u.segments)):
        seg = u.segments[seg_i]
        out_seg_prefix = seg.segid[0:3].strip()

        num_chunks = int(len(seg.residues) / args.chunk_size)
        if len(seg.residues) % args.chunk_size != 0:
            num_chunks += 1

        if num_chunks == 1:
            out_segid = seg.segid
            # We have to say delete=False here because the temporary file's variable may be GCed before we invoke VMD below,
            # which will cause the file to be deleted
            tmp = tempfile.NamedTemporaryFile(prefix=f'{out_segid}_', suffix='.pdb', delete=False)
            seg.atoms.write(tmp.name)
            print(f"Filename for segid #{seg_i} for {seg.segid} is {tmp.name}", file=sys.stderr)
            segs_to_process.append((out_segid, tmp.name))
            continue

        # Sanity check for fabulously large segment
        if num_chunks > len(out_seg_suffixes):
            print(f"You need {num_chunks} chunks for segment {seg.segid} but I am only smart enough for {len(out_seg_suffixes)}", file=sys.stderr)
            exit(1)

        for chunk_i in range(num_chunks):
            out_segid = f"{out_seg_prefix}{out_seg_suffixes[chunk_i]}"

            length = args.chunk_size
            offset = chunk_i * length
            if (offset + length) >= len(seg.residues):
                length = len(seg.residues) - offset
            
            # print(f"{seg.segid} (chunk #{chunk_i}) -> {out_segid}: residue offset = {offset}, length = {length}", file=sys.stderr)            
            tmp = tempfile.NamedTemporaryFile(prefix=f'{out_segid}_', suffix='.pdb', delete=False)
            print(f"Filename for segid #{seg_i} chunk #{chunk_i} for {seg.segid} is {tmp.name}", file=sys.stderr)

            # Apparently, we can't set the segid of an AtomGroup, but we can of a segment?
            # This does not work: seg.residues[offset:offset+length].atoms.segid = out_segid
            old_segid = seg.segid
            seg.segid = out_segid
            seg.residues[offset:offset+length].atoms.write(tmp.name)
            seg.segid = old_segid
            
            segs_to_process.append((out_segid, tmp.name))


    # Convert a newline-separated list of files to psfgen commands
    toppar_cmds = '\n'.join([f"topology {x.strip()}" for x in open(args.toppar_txt).readlines()])

    # Start generating that script
    psfgen_script = f"""package require psfgen
{toppar_cmds}
"""
    for segid, fname in segs_to_process:
        s = f"""
segment {segid} {{
    pdb {fname}
}}
coordpdb {fname} {segid}
"""
        psfgen_script += s

    # We let user specify output prefix
    psfgen_script += f"""
writepsf {args.out_prefix}.psf
writepdb {args.out_prefix}.pdb
"""

    if args.show_psfgen_script:
        print('psfgen script:', file=sys.stderr)
        print(psfgen_script)

    # Run vmd -dispdev text -eofexit
    print('Running VMD to do the actual PSF/PDB generation. Fasten your seatbelt...', file=sys.stderr)
    c = subprocess.run([args.vmdbin, '-dispdev', 'text', '-eofexit'], input=bytes(psfgen_script, 'utf-8'), capture_output=True)
    print(f'Done running VMD. I hope it worked! If any doubt, run me again with --show-vmd-output.', file=sys.stderr)
    
    if args.show_vmd_output:
        print('VMD stdout:\n===========', file=sys.stderr)
        print(str(c.stdout, 'utf-8'))
        print('\nVMD stderr:\n===========', file=sys.stderr)
        print(str(c.stderr, 'utf-8'))

    # Delete all those temporary files
    if args.keep_temp_pdbs is False:
        for segid, fname in segs_to_process:
            remove(fname)
        print('Deleted those temporary PDBs.', file=sys.stderr)
    else:
        print('Did not delete all those temporary PDBs.', file=sys.stderr)

if __name__ == "__main__":
    main()
