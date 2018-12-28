#!/usr/bin/env python
#
# Generates a blueprint file to be used as input to the Rosetta 3 "remodel" tool,
# so that it can generate missing loops, most likely with an ab initio protocol.

import argparse
import sys
import os
from shutil import copy2
import MDAnalysis as mda
from Bio import SeqIO

class PDBEmitter(object):
    """Stateful emitter for PDB (Protein Data Bank) formatted files.

    Takes MDAnalysis atoms. Automatically renumbers ATOMs but does not change residue IDs.
    This is intended to assist in merging two structures into one PDB file.

    Example:
        pdb = PDBEmitter(sys.stdout)
        u1 = MDAnalysis.Universe('foo.pdb')
        u2 = MDAnalysis.Universe('bar.pdb')
        pdb.emit_atom(u2.atoms[13])
        pdb.emit_atom(u1.atoms[42])
    """

    def __init__(self, out):
        self.out = out
        self.curr_atom_id = 1
        self.last_emitted_resid = 0
        self.last_encountered_resid = None

    def emit_atom(self, atom, occupancy=0.0, resid=None, renumber_resid=False):
        """Emits one ATOM record.

        Returns emitted resid."""
        if resid is None:
            if renumber_resid is True:
                resid = self.last_emitted_resid
                if atom.residue.resid != self.last_encountered_resid:
                    resid += 1
            else:
                resid = atom.residue.resid
        # Format string from: http://cupnet.net/pdb-format/
        self.out.write(
            '{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format(
            "ATOM", self.curr_atom_id, atom.name, atom.altLoc, atom.residue.resname, atom.residue.segid, resid,
                  ' ', atom.position[0], atom.position[1], atom.position[2],
                  occupancy, 0.0, atom.type, '  ')) # occupancy, temp factor, type, charge
        self.curr_atom_id += 1
        self.last_emitted_resid = resid
        self.last_encountered_resid = atom.residue.resid
        return resid


def main():
    ap = argparse.ArgumentParser(description='Generate Rosetta remodel blueprint file for ab initio loop modeling')
    ap.add_argument('pdb', help='PDB file of your protein with missing loops')
    ap.add_argument('segid', help='Which segid you want. Sorry, only one is allowed')
    ap.add_argument('fasta', help='FASTA sequence file of the specified segid')
    ap.add_argument('-o', '--output-dirname', help='Put output in this directory (current directory is not a good idea)')
    ap.add_argument('-l', '--loops', help='Care only about loops starting at these comma-separated resids (default: all)')
    ap.add_argument('-m', '--max-loop-length', type=int, default=50, help='Maximum length of loop to fill')
    args = ap.parse_args()

    # only_these_loops = [int(x) for x in args.loops.split(',')] if args.loops is not None else None

    u = mda.Universe(args.pdb)
    protein = u.select_atoms('protein and segid %s' % args.segid)
    sequence = SeqIO.read(args.fasta, 'fasta').seq

    output_dirname = os.path.abspath(args.output_dirname or 'remodel_%s_%s' % (args.pdb, args.segid))
    if not os.path.exists(output_dirname):
        os.mkdir(output_dirname)
    copy2(args.pdb, output_dirname)
    os.chdir(output_dirname)

    # Renumber PDB for Rosetta, because it will silently generate garbage unless resids start from
    # 1 and increase without gaps
    renumbered_protein_filename = 'renumbered-protein.pdb'
    with open(renumbered_protein_filename, 'w') as f:
        e = PDBEmitter(f)
        real_to_rosetta_resid = dict()
        for a in protein.atoms:
            rosetta_resid = e.emit_atom(a, renumber_resid=True)
            real_to_rosetta_resid[a.residue.resid] = rosetta_resid

        with open('real-to-rosetta-resid.map', 'w') as f2:
            for real_resid in sorted(real_to_rosetta_resid.keys()):
                f2.write('%d %d\n' % (real_resid, real_to_rosetta_resid[real_resid]))

    # First determine where all the loops are
    _, all_loop_starts, all_loop_lengths = gen_blueprint(protein, sequence)

    print >>sys.stderr, '##### Going back and generating input files for Rosetta remodel #####'

    # Now write individual blueprint and conf files for each loop
    # We do each loop separately so we can parallelize. It will be a little annoying to combine all the loops
    # back again but we can probably write a script to do it.
    for i in range(len(all_loop_starts)):
        loop_start, loop_length = all_loop_starts[i], all_loop_lengths[i]

        if loop_length <= args.max_loop_length:
            os.chdir(output_dirname)
            blueprint_entries, loop_starts, loop_lengths = gen_blueprint(protein, sequence, [loop_start,])
            assert loop_start == loop_starts[0], 'Bug: confused about loop start'
            assert loop_length == loop_lengths[0], 'Bug: confused about loop length'
            this_output_dirname = 'loop%d_%d' % (loop_start, loop_length)
            if not os.path.exists(this_output_dirname):
                os.mkdir(this_output_dirname)
            os.chdir(this_output_dirname)
            prefix = 'remodel'
            with open('%s.blueprint' % prefix, 'w') as f:
                print >>f, '\n'.join(blueprint_entries)
            with open('%s.conf' % prefix, 'w') as f:
                print >>f, """-in:file:s ../%s
-ignore_zero_occupancy false
-remodel:blueprint %s.blueprint
-run:chain %s
-remodel:num_trajectory 1
-nstruct 500
-out:level 100
-out:path:all .
-out:file:scorefile %s.sc
""" % (renumbered_protein_filename, prefix, args.segid, prefix)


# Generates blueprint file which tells Rosetta Remodel what to do with each residue in the protein
# (e.g. nothing, or generate a conformation)
def gen_blueprint(protein, sequence, only_these_loops=None):
    last_res = None
    rosetta_resid = 1
    blueprint_entries = []
    loop_starts, loop_lengths = [], []

    # Look for residues missing in the structure
    for res in protein.residues:
        res_1let = mda.lib.util.convert_aa_code(res.resname)
        res_seq = sequence[res.resid-1]
        # Ensure that the supplied sequence matches that in the PDB
        # Don't forget, resid counts starting at 1
        if res_1let != res_seq:
            print >>sys.stderr, 'Error: res_1let %s does not match sequence %s at %d. Is this the right FASTA? Is the PDB numbered properly?' % \
                (res_1let, res_seq, res.resid)
            exit(1)

        if last_res is None:
            last_res = res

        discon_length = res.resid - last_res.resid - 1
        if discon_length > 0 and (only_these_loops is None or last_res.resid in only_these_loops):
            print >>sys.stderr, '# Discontinuity between residues %s:%s%d and %s:%s%d (length %d)' \
                   % (last_res.segid, last_res.resname, last_res.resid % 10000,
                      res.segid, last_res.resname, res.resid % 10000, discon_length)

            # We need to include the residues on either end of the loop, for which we do know
            # the coordinates, in the loop modeling, so as to allow sufficient flexibility
            lastres_1let = mda.lib.util.convert_aa_code(last_res.resname)
            # The last entry in the blueprint must be replaced
            blueprint_entries.pop()
            # Tell Rosetta remodel to pick the correct amino acid ("PIKAA" command) for this _L_oop
            blueprint_entries.append('%d %s L PIKAA %s' % (rosetta_resid-1, lastres_1let, lastres_1let))
            for i in range(last_res.resid+1, res.resid):
                blueprint_entries.append('0 X L PIKAA %s' % sequence[i-1])
            blueprint_entries.append('%d %s L PIKAA %s' % (rosetta_resid, res_1let, res_1let))

            # blueprint_entries.append('%d %s L PIKAA %s' % (last_res.resid, lastres_1let, lastres_1let))
            # for i in range(last_res.resid+1, res.resid):
            #     blueprint_entries.append('0 X L PIKAA %s' % sequence[i-1])
            # blueprint_entries.append('%d %s L PIKAA %s' % (res.resid, res_1let, res_1let))
            loop_starts.append(last_res.resid)
            loop_lengths.append(discon_length)
        else:
            # Tell Rosetta remodel to do nothing ("." command)
            blueprint_entries.append('%d %s .' % (rosetta_resid, res_1let))
            # blueprint_entries.append('%d %s .' % (res.resid, res_1let))

        last_res = res
        rosetta_resid += 1

    return blueprint_entries, loop_starts, loop_lengths


if __name__ == '__main__':
    main()
