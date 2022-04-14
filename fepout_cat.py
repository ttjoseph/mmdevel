#!/usr/bin/env python3
#
# Concatenates a bunch of NAMD fepout files, assuming they all have the same lambdas.
# Why not use alchemlyb for this? Because alchemlyb insists on a higher level of self-
# consistency in the fepouts. Here, all fepouts *must* be of the same window. The idea
# is to preprocess a mess of fepouts, some of which won't have any lambda values, so
# that alchemlyb will accept them. I don't want to put this funcitonality in alchemlyb
# because it is not intuitive and will lead to user error and corrupted (possibly
# published) results.
#
# Because this script does not actually do very much, we can take a dumb approach to
# parsing the fepouts. We will pass through all the lines except for those
#
# Tom Joseph, University of Pennsylvania
import re
import argparse
import sys
import fileinput

def main():
    ap = argparse.ArgumentParser(description='Concatenates fepout files from a single lambda window')
    ap.add_argument('--equil-steps', '-e', type=int, help='Number of equilibration steps to assume')
    ap.add_argument('filenames', nargs='+', help='.fepout file names')
    args = ap.parse_args()

    markers_re = {'new_fep_window': re.compile(r'#NEW FEP WINDOW: LAMBDA SET TO ([\d.]+) LAMBDA2 ([\d.]+)'),
        # 'new_fep_window_idws': re.compile(r'#NEW FEP WINDOW: LAMBDA SET TO ([\d.]+) LAMBDA2 ([\d.]+) LAMBDA_IDWS ([\d.]+)'),
        'equil_end': re.compile(r'#\d+ STEPS OF EQUILIBRATION AT LAMBDA [\d.]+ COMPLETED'),
        'prod_start': re.compile(r'#STARTING COLLECTION OF ENSEMBLE AVERAGE'),
        'prod_end': re.compile(r'#Free energy change for lambda window \[ ([\d.]+) ([\d.]+) \] is ([-\deE.]+) ; net change until now is [-\deE.]+')}
    markers_line = dict.fromkeys(markers_re.keys())
    comment_lines = set()
    prev_step_count = None # Used to guess fepout output frequency
    output_frequency = None
    num_data_lines = 0
    equil_done = False

    # The fileinput module succinctly iterates through a bunch of files as if they were one,
    # and transparently opening compressed files
    with fileinput.input(files=sorted(args.filenames), openhook=fileinput.hook_compressed) as f:
        for line in f:
            # Check if we've seen any of the expected regexes. If so, save the line so we can
            # ensure the the rest of said line is consistent.
            matched = False
            for key, regex in markers_re.items():
                m = regex.match(line)
                if m:
                    matched = True
                    if key == 'new_fep_window':
                        lambda1, lambda2 = m.group(1), m.group(2)
                    if markers_line[key] == None:
                        markers_line[key] = line
                        should_print_marker = True
                        # Wait till end of all the files to print prod_end marker
                        if key == 'prod_end':
                            should_print_marker = False
                        # Don't print equil_end/prod_start markers if user specified equil-steps
                        if args.equil_steps is not None and key in ('equil_end', 'prod_start'):
                            should_print_marker = False
                        if should_print_marker:
                            print(line, end='')
                        break
                    else:
                        # We've seen this type of marker before.
                        # TODO: If this marker has lambda values, check that they haven't changed.
                        pass

            # TODO:
            # Are we going to manually insert the equil_end and prod_start markers? If not, just pass through
            # everything but prod_end
            if output_frequency != None:
                num_steps = num_data_lines * output_frequency
                # If user specified the number of equilibraitions steps we *should* have, enforce that
                # by outputting a marker line at the right time
                if not equil_done and args.equil_steps is not None and num_steps >= args.equil_steps:
                    print(f'{args.equil_steps} STEPS OF EQUILIBRATION AT LAMBDA {lambda1} COMPLETED')
                    print('#STARTING COLLECTION OF ENSEMBLE AVERAGE')
                    equil_done = True


            # Perhaps this is a comment line - filter out duplicates
            if not matched and line.startswith('#') and line not in comment_lines:
                print(line, end='')
                comment_lines.add(line)
            elif not matched:
                # Regular old data line.
                # Keep track of how many of these we've seen so we can decide what to emit when (above).
                # Guess the output frequency if we haven't already.                
                if line.startswith('FepE'):
                    num_data_lines += 1
                    if output_frequency == None:
                        tokens = line.split()
                        step_count = int(tokens[1])
                        if prev_step_count == None:
                            prev_step_count = step_count
                        else:
                            putative_output_frequency = step_count - prev_step_count
                            if putative_output_frequency > 0:
                                output_frequency = putative_output_frequency
                                print(f'{fileinput.filename()}: Guessed output frequency to be', output_frequency, file=sys.stderr)


                print(line, end='')
        print(markers_line['prod_end'], end='')


if __name__ == '__main__':
    main()