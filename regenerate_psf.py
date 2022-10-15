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
from os.path import isfile, getmtime
import tempfile
import subprocess
from libttj import ShadyPDB, ShadyPSF, renumber_psf_pdb

def main():
    ap = argparse.ArgumentParser(description='Remakes a PSF/PDB from existing PDB. Helpful if you want to avoid the CHARMM-GUI rigamarole')
    ap.add_argument('pdb', help='PDB file to use')
    ap.add_argument('toppar_txt', help='File containing newline-separated list of psfgen-acceptable topology files')
    ap.add_argument('out_prefix', help='Prefix for output PSF and PDB files')
    ap.add_argument('--disulfide-bonds', help='Specify disulfide bonds. Ex: PROA:1-PROA:2,PROA:4-PROA:5')
    ap.add_argument('--patch', help='Specify arbitrary single-residue patch. Ex: PROA:1=SP2,PROB:42=SP1')
    ap.add_argument('--mutate', help='Specify residue mutations. Ex: PROA:1=FOO,PROB:42=BAR')
    ap.add_argument('--hmassrepart', action='store_true', help='Tell psfgen to do hydrogen mass repartitioning')
    ap.add_argument('--chunk-size', type=int, default=9999, help='Maximum number of residues in a segment that you think psfgen or whatever can handle')
    ap.add_argument('--show-psfgen-script', action='store_true', help='Show the psfgen script that will be passed to VMD')
    ap.add_argument('--vmdbin', default='vmd', help='VMD binary')
    ap.add_argument('--keep-temp-pdbs', action='store_true', help='Keep the temporary PDBs made for psfgen so you can look at them')
    args = ap.parse_args()

    # psfgen needs to have input PDBs chopped up into segments and regurgitated into it, much like a helpless baby bird.
    # I can't figure out how to get VMD to deal with a number of residues > 9999 when we start with a PDB so we will do the segment
    # splitting here, in Python. (This often is necessary when there is a bunch of water in a simulation system.)
    out_seg_suffixes = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    u = ShadyPDB(args.pdb)
    segs_to_process = []
    seg_i = 0

    # Iterate over tuples that associate each residue index (NOT resid) to each segid.
    # We use residue *index* rather than resid because the PDB resid field is too small
    # to accommodate a lot of residues.
    for segid, resindexes in u.segid_to_resindex.items():
        out_seg_prefix = segid[0:3].strip()

        num_chunks = int(len(resindexes) / args.chunk_size)
        if len(resindexes) % args.chunk_size != 0:
            num_chunks += 1

        if num_chunks == 1:
            out_segid = segid
            # We have to say delete=False here because the temporary file's variable may be GCed before we invoke VMD below,
            # which will cause the file to be deleted
            tmp = tempfile.NamedTemporaryFile(prefix=f'{out_segid}_', suffix='.pdb', mode='wt', delete=False)
            # Replacement for seg.atoms.write(tmp.name)
            chunk_atomindexes = []
            for resindex in resindexes:
                chunk_atomindexes.extend(u.resindex_to_atomindex[resindex])
            print(f"Filename for segid #{seg_i} for {segid} ({len(chunk_atomindexes)} atoms) is {tmp.name}", file=sys.stderr)
            # Preserve residue IDs so patch commands act as expected by user
            u.renumber_subset(chunk_atomindexes, renumber_residues=False)
            u.write_to_pdb(tmp, chunk_atomindexes)
            segs_to_process.append((out_segid, tmp.name))
            seg_i += 1
            continue

        # Sanity check for fabulously large segment where we've run out of segment name suffix characters
        if num_chunks > len(out_seg_suffixes):
            print(f"Error: You need {num_chunks} chunks for segment {segid} but I am only smart enough for {len(out_seg_suffixes)}", file=sys.stderr)
            exit(1)

        # We would only reach here if num_chunks > 1 (note continue statement above)
        for chunk_i in range(num_chunks):
            out_segid = f"{out_seg_prefix}{out_seg_suffixes[chunk_i]}"

            length = args.chunk_size
            offset = chunk_i * length
            if (offset + length) >= len(resindexes):
                length = len(resindexes) - offset
            
            # print(f"{seg.segid} (chunk #{chunk_i}) -> {out_segid}: residue offset = {offset}, length = {length}", file=sys.stderr)
            print(f"offset: {offset}; length: {length}; len(resindexes): {len(resindexes)}",
                file=sys.stderr)
            tmp = tempfile.NamedTemporaryFile(prefix=f'{out_segid}_', suffix='.pdb', mode='wt', delete=False)

            # Get all the atom indexes for this chunk then process them
            # seg.residues[offset:offset+length].atoms.write(tmp.name)
            # Write this subset of atoms with a custom segid
            # offset is relative to the first residue in this segment, where offset 0 is the first residue
            chunk_atomindexes = []
            for resindex in resindexes[offset:offset+length]:
                chunk_atomindexes.extend(u.resindex_to_atomindex[resindex])
            
            print(f"Filename for segid #{seg_i} chunk #{chunk_i} for {segid} ({len(chunk_atomindexes)} atoms) is {tmp.name}", file=sys.stderr)
            u.renumber_subset(chunk_atomindexes)
            u.set_segid(chunk_atomindexes, out_segid)
            u.write_to_pdb(tmp, chunk_atomindexes)
            
            segs_to_process.append((out_segid, tmp.name))
        seg_i += 1

    # Convert a newline-separated list of files to psfgen commands
    toppar_lines = open(args.toppar_txt).readlines()
    toppar_lines = [l for l in (line.strip() for line in toppar_lines) if l]
    toppar_cmds = '\n'.join([f"topology {x.strip()}" for x in toppar_lines])

    out_psf_filename, out_pdb_filename = f'{args.out_prefix}.psf', f'{args.out_prefix}.pdb'
    
    # Start generating that script
    psfgen_script = f"""package require psfgen
{toppar_cmds}
"""
    for segid, fname in segs_to_process:
        # Grab segid-specific mutate commands
        # Ex: PROA:3=ASP -> mutate 3 ASP
        # But only if segd == PROA
        mutate_str = ''
        if args.mutate != None:
            for mutation in args.mutate.split(','):
                mutate_from, mutate_to = mutation.split('=')
                mutate_segid, mutate_resid = mutate_from.split(':')
                if segid == mutate_segid:
                    mutate_str += f'mutate {mutate_resid} {mutate_to}\n'

        psfgen_script +=f"""
segment {segid} {{
    pdb {fname}
    {mutate_str}
}}
"""
    # Add disulfide bonds
    if args.disulfide_bonds is not None:
        for disu in args.disulfide_bonds.split(','):
            res1, res2 = disu.split('-')
            psfgen_script += f"patch DISU {res1} {res2}\n"

    # Add any arbitrary single-residue patches
    if args.patch is not None:
        for patch_str in args.patch.split(','):
            res, patch = patch_str.split('=')
            psfgen_script += f"patch {patch} {res}\n"

    # Specify coordinates after the patches for disulfide bonds, according to psfgen user guide example
    for segid, fname in segs_to_process:
        psfgen_script += f"coordpdb {fname} {segid}\n"

    if args.patch is not None or args.mutate is not None:
        psfgen_script += "# Guess missing coordinates resulting from patches/mutations\nguesscoord\n"
        psfgen_script += "# And we need to regenerate angle and dihedral parameters\nregenerate angles dihedrals"

    if args.hmassrepart:
        psfgen_script += "hmassrepart 1\n"

    # We let user specify output prefix
    psfgen_script += f"""
writepsf {out_psf_filename}
writepdb {out_pdb_filename}
"""

    if args.show_psfgen_script:
        print('psfgen script:', file=sys.stderr)
        print(psfgen_script)

    # Remember modified time of PSF and PDB files if they exist, because we're going to overwrite them.
    # If there was an error, we want to know that 
    psf_mtime = getmtime(out_psf_filename) if isfile(out_psf_filename) else None
    pdb_mtime = getmtime(out_pdb_filename) if isfile(out_pdb_filename) else None        
    
    # Run vmd -dispdev text -eofexit
    print('Running VMD to do the actual PSF/PDB generation. Fasten your seatbelt...', file=sys.stderr)
    c = subprocess.run([args.vmdbin, '-dispdev', 'text', '-eofexit'], input=bytes(psfgen_script, 'utf-8'), capture_output=True)
    print(f'Done running VMD. I hope it worked!', file=sys.stderr)
    
    # Check if VMD actually made the PSF/PDB we asked for
    should_renumber = True
    if isfile(out_psf_filename) is False:
        print(f'Error: The PSF file that was supposed to be generated does not appear to exist.')
        should_renumber = False
    elif getmtime(out_psf_filename) == psf_mtime:
        print(f'Error: Does not look like the PSF file was actually regenerated by VMD.')
        should_renumber = False

    if isfile(out_pdb_filename) is False:
        print(f'Error: The PDB file that was supposed to be generated does not appear to exist.')
        should_renumber = False
    elif getmtime(out_pdb_filename) == pdb_mtime:
        print(f'Error: Does not look like the PDB file was actually regenerated by VMD.')
        should_renumber = False
    
    vmd_output_file = f'{args.out_prefix}.vmd.log'
    with open(vmd_output_file, 'w') as f:
        print('VMD stdout:\n===========', file=f)
        print(str(c.stdout, 'utf-8'), file=f)
        print('\nVMD stderr:\n===========', file=f)
        print(str(c.stderr, 'utf-8'), file=f)
        print(f'Wrote VMD output to {vmd_output_file}.')

    # Delete all those temporary files
    if args.keep_temp_pdbs is False:
        for segid, fname in segs_to_process:
            remove(fname)
        print('Deleted those temporary PDBs.', file=sys.stderr)
    else:
        print('Did not delete all those temporary PDBs.', file=sys.stderr)

    if should_renumber:
        print(f'Renumbering the solute residues in the {out_psf_filename} and {out_pdb_filename} to match {args.pdb}...', file=sys.stderr)
        psf = ShadyPSF(out_psf_filename)
        pdb = ShadyPDB(out_pdb_filename)
        pdb_template = ShadyPDB(args.pdb)
        psf, pdb = renumber_psf_pdb(psf, pdb, pdb_template)
        print(f'Writing {out_psf_filename}...', file=sys.stderr)
        psf.write_to_psf(out_psf_filename)
        print(f'Writing {out_pdb_filename}...', file=sys.stderr)
        pdb.write_to_pdb(out_pdb_filename)
        print(f'All done. Enjoy!', file=sys.stderr)

if __name__ == "__main__":
    main()
