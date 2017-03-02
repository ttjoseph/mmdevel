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
import cStringIO as StringIO
import yaml
from natsort import natsorted
import intervaltree

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

def trim_timesteps_from_block(block, numsteps=0):
    """Trims numsteps timesteps from the start of the provided FepEnergy block.

    Note this is insensitive to the the interval of the steps. The same amount of
    simulation time is trimmed regardless. So if you wanted 100000 steps trimmed but
    had an output interval of 20000 (for some reason), only five FepEnergy values
    would be trimmed."""

    in_buf = StringIO.StringIO(block)
    out_buf = StringIO.StringIO()
    first_step = None
    for line in in_buf.readlines():
        if line.startswith('FepEnergy:'):
            step = int(line.split()[1])
            if first_step is None: first_step = step
            if (step - first_step) >= numsteps:
                out_buf.write(line)
        else:
            out_buf.write(line)

    new_block = out_buf.getvalue()
    out_buf.close()
    return new_block


def concat_block_prod(blocks, first_timestep=0):
    """Returns a string containing combined production parts of the given blocks, including whatever header
    is necessary for VMD ParseFEP to be happy."""

    # Load the raw data content of the block
    blocks_data = [get_block_from_file(b.data['fname'], b.data['prod_start_offset'], b.data['prod_end_offset']) for b in blocks]
    out_buf = StringIO.StringIO()
    headers, fepenergy_data, footers = [], [], []
    delta = None

    # Get the header from each block and ensure they are all for the same lambda range
    for block_data in blocks_data:
        in_buf = StringIO.StringIO(block_data)
        header, header_done, footer = [], False, []
        fepenergy = []
        for line in in_buf.readlines():
            # Parse headers and footers, which all start with '#'
            if line.startswith('#'):
                if header_done is False:
                    header.append(line)
                else:
                    footer.append(line)
            elif line.startswith(('FepEnergy:')):
                fepenergy.append(line)
                # If we see a FepEnergy line, we are by definition past the header
                header_done = True

        headers.append(header)
        fepenergy_data.append(fepenergy)
        footers.append(footer)

    # Ensure the lambda ranges are the same, by induction. TODO this doesn't work because it's not always
    # in the third line.
    # Also guess the FepEnergy delta
    # for block_i in range(1, len(headers)):
    #     #NEW FEP WINDOW: LAMBDA SET TO %f LAMBDA2 %f
    #     if headers[block_i][2] != headers[block_i-1][2]:
    #         sys.stderr.write('concat_block_prod: The lambdas in these blocks do not match:\n')
    #         sys.stderr.write(headers[block_i-1][2])
    #         sys.stderr.write('\n')
    #         sys.stderr.write(headers[block_i][2])
    #         sys.stderr.write('\n')
    #         return None

    # Keep the header from the first block
    for line in headers[0]:
        out_buf.write(line)

    # Concatenate the FepEnergy lines, renumbering as appropriate
    fepenergy_timestep = first_timestep
    for block_i in range(len(fepenergy_data)):
        # FepEnergy: 50 ...
        # FepEnergy: 100 ...
        this_delta = int(fepenergy_data[block_i][1].split()[1]) - int(fepenergy_data[block_i][0].split()[1])

        for line in fepenergy_data[block_i]:
            fepenergy_timestep += this_delta
            # FepEnergy: %6d %14.4f %14.4f...
            out_buf.write('FepEnergy: %6d' % fepenergy_timestep)
            fepenergy_timestep += this_delta
            tokens = line.strip().split()
            for i in range(2, len(tokens)):
                out_buf.write(' %14.4f' % float(tokens[i]))
            out_buf.write('\n')

    # Preserve the footer from the last block
    for line in footers[-1]:
        out_buf.write(line)
        out_buf.write('\n')

    # Return the assembled concatenated block
    new_block = out_buf.getvalue()
    out_buf.close()
    return new_block


def load_spec_file(specfile):
    spec = intervaltree.IntervalTree()
    spec_fepout_files = set()
    spec_yaml = yaml.load(open(specfile))

    # Parse the spec file and save it in an interval tree
    for lambda_range in spec_yaml:
        fepouts = spec_yaml[lambda_range]
        (spec_b, spec_e) = sorted(float(x) for x in lambda_range.strip().split())
        # If we haven't seen this interval before, store it
        if not spec.search(spec_b, spec_e, strict=True):
            # ...But we may overlap with a previous interval, which is an error
            if spec.search(spec_b, spec_e):
                sys.stderr.write(
                    '(%f, %f) overlaps with something else in the spec. I am too lazy to tell you what\n' % (
                    spec_b, spec_e))
                return None, None
            spec[spec_b:spec_e] = []
        current = list(spec[spec_b:spec_e])[0].data
        current.append(fepouts)
        spec[spec_b:spec_e] = current
        spec_fepout_files.update(set(fepouts))
        # sys.stderr.write('    %s: lambda %f to %f\n' % (fepout, b, e))

    return spec, spec_fepout_files

def main():
    ap = argparse.ArgumentParser(description='Stitch together NAMD fepout files, discarding incomplete lambda windows, and dump the result to stdout')
    ap.add_argument('--spec', help='FEP specification that dictates which lambda ranges should be used from which fepout files')
    ap.add_argument('--ignore-steps', type=int, default=0, help='Ignore this many steps at the start of each production window')
    ap.add_argument('fepout_file', nargs='+', help='.fepout files, in order of lambdas you want')
    args = ap.parse_args()

    # The general strategy is to record offsets within files delineating blocks we want to keep.
    # This way we can avoid keeping big blocks of data around in memory.
    last_delta = 0

    spec, spec_fepout_files = None, set()

    # If user provides a yaml file that describes which lambda ranges can come from which files,
    # go through each block that we inhaled and keep only those which meet the criteria.
    # Bonus: Maybe check to ensure that we have covered the entire transformation.
    if args.spec:
        if not isfile(args.spec):
            sys.stderr.write('Cannot find your FEP spec file %s.' % args.spec)
            return 1

        spec, spec_fepout_files = load_spec_file(args.spec)
        
    # If user provided a spec, there's no need to load fepout files that aren't in the spec.
    # So get rid of all fepout files in our list that we can't need.
    all_fepout_files = set(args.fepout_file)
    if len(spec_fepout_files) > 0:
        all_fepout_files = all_fepout_files & spec_fepout_files

    # Load all blocks from all fepout files and store them in an interval tree
    all_blocks = intervaltree.IntervalTree()
    all_lambda_ranges = set()
    for fname in all_fepout_files:
        # Get a list of blocks in this fepout file
        blocks = scan_fepout_file(fname)

        for key in blocks:
            (block_b, block_e) = sorted(float(l) for l in key.split('_'))
            all_lambda_ranges.add((block_b, block_e))

            # Look through all lambda ranges in spec to see whether this block is allowed
            if spec:
                spec_files = spec.search(block_b, block_e)
                if len(spec_files) == 0: # Not in spec, keep on truckin'
                    continue
                elif len(spec_files) > 1:
                    sys.stderr.write('(%f, %f) spans more than one spec file entry. I am confused.\n' % (block_b, block_e))
                    return 1
                else: # We found exactly one result so keep going
                    sys.stderr.write('Keeping key (%f, %f)\n' % (block_b, block_e)) # TODO: More detail
                    pass

            # Dump this block into the tree. Duplicates are OK because they'll be returned with intervaltree.search()
            all_blocks[block_b:block_e] = blocks[key]

            # We use the last delta value encountered to decide whether lambdas increase or decrease.
            # Naturally this does not protect against weird pathological cases where the user provides
            # both increasing and decreasing lambda windows. I guess we should error out in that case,
            # but such intelligence is not yet implemented. Could easily be done by comparing the sign
            # of last_delta and offsets[key]['delta']; if they are different, it's bad news bears.
            last_delta = blocks[key]['delta']

    # Spit out a single header block from any .fepout file
    for block in all_blocks:
        sys.stdout.write(get_block_from_file(block.data['fname'], 0, block.data['header_end_offset']))
        break

    total_energy_change = 0.0

    # TODO:
    for (block_b, block_e) in all_lambda_ranges:
        # Get all blocks with this lambda range
        blocks = all_blocks.search(block_b, block_e)
        print >>sys.stderr, concat_block_prod(blocks)
        sys.stderr.write('%s from %s: Equil: %d bytes; Prod %d bytes\n' % (key,
            block.data['fname'],
            block.data['equil_end_offset'] - block.data['equil_start_offset'],
            block.data['prod_end_offset'] - block.data['prod_start_offset']))
        sys.stdout.write(get_block_from_file(block.data['fname'],
            block.data['equil_start_offset'],
            block.data['equil_end_offset']))
        prod_block = get_block_from_file(block.data['fname'],
            block.data['prod_start_offset'],
            block.data['prod_end_offset'])

        prod_block = trim_timesteps_from_block(prod_block, args.ignore_steps)
        sys.stdout.write(prod_block)

        total_energy_change += block.data['energy_change']

    sys.stderr.write('Kept %d lambda windows.\n' % len(all_blocks))
    sys.stderr.write('Total energy change: %f\n' % total_energy_change)
    return 0

if __name__ == '__main__':
    exit(int(main() or 0))
