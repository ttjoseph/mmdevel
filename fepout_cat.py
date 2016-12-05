#!/usr/bin/env python
#
# This script is intended to combine multiple .fepout files from a single NAMD FEP calculation
# that was prematurely terminated (perhaps due to cluster queue time limits) and then restarted.
#
# We make the following assumptions in this script:
# - The beginning of each fepout file is not corrupted
#
# Input: Multiple .fepout files, in lambda order
# Output: Single .fepout, sent to stdout.

import re
import argparse
import sys
from os.path import isfile
import yaml
from natsort import natsorted

def readline_and_offset(f):
    offset = f.tell()
    return (f.readline(), offset)


def readline_and_offset_until(f, until_re):
    # Read a line. Does it match the regex? If so, return the line and its offset in the file.
    # If not, try the next line. If we hit EOF, also stop looking.
    while True:
        line, offset = readline_and_offset(f)
        if until_re.match(line) or line == '':
            return line, offset


def scan_fepout_file(fname):
    new_fep_window_re = re.compile(r'#NEW FEP WINDOW: LAMBDA SET TO ([\d.]+) LAMBDA2 ([\d.]+)')
    equil_end_re = re.compile(r'#\d+ STEPS OF EQUILIBRATION AT LAMBDA [\d.]+ COMPLETED')
    prod_start_re = re.compile(r'#STARTING COLLECTION OF ENSEMBLE AVERAGE')
    prod_end_re = re.compile(r'#Free energy change for lambda window \[ ([\d.]+) ([\d.]+) \] is ([-\d.]+) ; net change until now is [-\d.]+')

    f = open(fname)
    # Eat the start of file header, which should span two lines
    readline_and_offset(f)
    readline_and_offset(f)
    header_end_offset = f.tell()
    lambdas = {}
    while True:
        line, offset = readline_and_offset_until(f, new_fep_window_re)
        if line == '': return lambdas # Give up on EOF
        # We have started a new window so extract the current lambdas.
        # This is the equilibration stage.
        m = new_fep_window_re.match(line)
        l0, l1 = float(m.group(1)), float(m.group(2))
        key = '%f_%f' % (l0, l1)
        current_lambda = {'fname': fname, 'header_end_offset': header_end_offset, 'equil_start_offset': offset,
            'delta': l1 - l0}
        line, offset = readline_and_offset_until(f, equil_end_re)
        if line == '': return lambdas
        current_lambda['equil_end_offset'] = f.tell()
        
        # Now here is the production part.
        line, offset = readline_and_offset_until(f, prod_start_re)
        if line == '': return lambdas
        current_lambda['prod_start_offset'] = offset
        line, offset = readline_and_offset_until(f, prod_end_re)
        if line == '': return lambdas
        current_lambda['prod_end_offset'] = f.tell()
        # Save the free energy change for this window so we can do some calculations with it later
        current_lambda['energy_change'] = float(prod_end_re.match(line).group(3))
        lambdas[key] = current_lambda

def get_block_from_file(fname, start, end):
    with open(fname) as f:
        f.seek(start)
        return f.read(end - start)

def main():
    ap = argparse.ArgumentParser(description='Stitch together NAMD fepout files, discarding incomplete lambda windows, and dump the result to stdout')
    ap.add_argument('--spec', help='FEP specification that dictates which lambda ranges should be used from which fepout files')
    ap.add_argument('fepout_file', nargs='+', help='.fepout files, in order of lambdas you want')
    args = ap.parse_args()

    # The general strategy is to record offsets within files delineating blocks we want to keep.
    # This way we can avoid keeping big blocks of data around in memory.
    all_offsets = {}
    last_delta = 0
    spec = None

    # TODO:
    # If user provides a yaml file that describes which lambda ranges can come from which files,
    # go through each block that we inhaled and keep only those which meet the criteria.
    # Bonus: Maybe check to ensure that we have covered the entire transformation.
    if args.spec:
        if not isfile(args.spec):
            sys.stderr.write('Cannot find FEP your spec file %s.' % args.spec)
            return 1
        
        spec_yaml = yaml.load(open(args.spec))
        spec = {}
        sys.stderr.write('Here is the FEP specification you provided:\n')
        for fepout in spec_yaml:
            (b, e) = (float(l) for l in spec_yaml[fepout].split())
            spec[fepout] = (b, e)
            sys.stderr.write('    %s: lambda %f to %f\n' % (fepout, b, e)) 

    for fname in args.fepout_file:
        offsets = scan_fepout_file(fname)
        # Allow newer blocks to override older ones
        for key in offsets:
            # Look for a reason to not include this block.
            # First up: if this fepout file was not mentioned in the spec, toss the block
            if spec and fname not in spec:
                sys.stderr.write('Ignoring %s from %s because that fepout file was not in the spec\n' % (key, fname))
                continue
            # Next: this fepout file is mentioned in the spec, but this lambda window is disallowed
            if spec and fname in spec:
                (spec_b, spec_e) = spec[fname]
                # Extract the lambda window values from the block key
                (b, e) = (float(l) for l in key.split('_'))
                # This lambda window must be contained within the spec lambda range
                if b >= spec_b and e <= spec_e:
                    sys.stderr.write('Allowing %s: %f to %f because we want to keep windows in range %f to %f\n' % (fname, b, e,
                        spec_b, spec_e))
                else:
                    sys.stderr.write('Ignoring %s: %f to %f because that lambda window is not within the allowed range\n' % (fname, b, e))
                    continue

            all_offsets[key] = offsets[key]
            # We use the last delta value encountered to decide whether lambdas increase or decrease.
            # Naturally this does not protect against weird pathological cases where the user provides
            # both increasing and decreasing lambda windows. I guess we should error out in that case,
            # but such intelligence is not yet implemented. Could easily be done by comparing the sign
            # of last_delta and offsets[key]['delta']; if they are different, it's bad news bears.
            last_delta = offsets[key]['delta']

    # Spit out a single header block from any .fepout file
    for key in all_offsets:
        sys.stdout.write(get_block_from_file(all_offsets[key]['fname'], 0, all_offsets[key]['header_end_offset']))
        break

    total_energy_change = 0.0

    for key in natsorted(all_offsets.keys(), reverse=(last_delta<0)):
        block_info = all_offsets[key]
        sys.stderr.write('%s from %s: Equil: %d bytes; Prod %d bytes\n' % (key,
            block_info['fname'],
            block_info['equil_end_offset'] - block_info['equil_start_offset'],
            block_info['prod_end_offset'] - block_info['prod_start_offset']))
        sys.stdout.write(get_block_from_file(block_info['fname'],
            block_info['equil_start_offset'],
            block_info['equil_end_offset']))
        sys.stdout.write(get_block_from_file(block_info['fname'],
            block_info['prod_start_offset'],
            block_info['prod_end_offset']))

        total_energy_change += block_info['energy_change']

    sys.stderr.write('Kept %d lambda windows.\n' % len(all_offsets))
    sys.stderr.write('Total energy change: %f\n' % total_energy_change)
    return 0

if __name__ == '__main__':
    exit(int(main() or 0))