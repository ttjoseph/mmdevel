#!/usr/bin/env python3
import argparse
import sys
import re
from yaml import safe_load, tokens
import os.path
from glob import glob
from string import ascii_lowercase
from functools import reduce

# Given a description of the overall FEP calculation and a particular window's .fepout files,
# generate a new .namd file that continues the calculation for this window.
# Of note, the .fepout files may be truncated.

# Config file example (YAML format):
# prefix: dis5A
# idws: true
# lambdas: 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
# equil_steps_per_window: 200000
# prod_steps_per_window: 1000000


# We'll use this to determine whether the fepout files for this window have enough steps
def num_steps_in_fepout(fname):
    with open(fname) as f:
        min_step_num, max_step_num = float('inf'), -1
        a_few_steps = []
        for line in f:
            if not line.startswith('FepE'):
                continue
            tokens = line.split()
            step_num = int(tokens[1])
            # Keep track of maximum and minimum step numbers because we might not have started at zero
            if step_num > max_step_num:
                max_step_num = step_num
            if step_num < min_step_num:
                min_step_num = step_num
            # We'll calculate the difference between the first two STEP number values as the
            # alchemy output frequency
            if len(a_few_steps) < 2:
                a_few_steps.append(step_num)

    # Avoid fencepost error because the FepEnergy for the first timestep is not recorded
    total_steps = max_step_num - min_step_num + (a_few_steps[1] - a_few_steps[0])
    return total_steps


def lambdas_for_fepout(fname):
    with open(fname) as f:
        lambda1, lambda2, lambda_idws = None, None, None
        # Unfortunately we scan through the whole file to ensure there's only one lambda window in it
        for line in f:
            if line.startswith('#NEW FEP WINDOW'):
                # More than one lambda window in this file?
                # Not supported yet, because that's not how I do things
                if lambda1 is not None:
                    print(f"Error: More than one window in the same file {fname} - I don't know how to support this yet",
                        file=sys.stderr)
                    exit(1)
                # #NEW FEP WINDOW: LAMBDA SET TO 0.1 LAMBDA2 0.108333 LAMBDA_IDWS 0.091667
                tokens = line.split()
                assert(tokens[3] == 'LAMBDA')
                lambda1, lambda2 = float(tokens[6]), float(tokens[8])
                if 'LAMBDA_IDWS' in line:
                    lambda_idws = float(tokens[10])

        return lambda1, lambda2, lambda_idws

# From: https://stackoverflow.com/questions/48983939/convert-a-number-to-excel-s-base-26
def divmod_base26(n):
    a, b = divmod(n, 26)
    if b == 0:
        return a - 1, b + 26
    return a, b


def to_base26(num):
    chars = []
    while num > 0:
        num, d = divmod_base26(num)
        chars.append(ascii_lowercase[d - 1])
    return ''.join(reversed(chars))


def from_base26(chars):
    return reduce(lambda r, x: r * 26 + x + 1, map(ascii_lowercase.index, chars), 0)


# Returns True if any pathname in list filenames is a file
def any_isfile(filenames):
    for fname in filenames:
        if os.path.isfile(fname):
            return True
    return False

# Strategy for generating the 'next' name:
# First part is common prefix
# Second part is a zero-padded 3-digit number
# Third part one of 'a', 'b', ..., 'aa', 'ab', ...according to the base-26 column numbering scheme
def generate_next_prefix(config, fnames):
    middles = [x.replace(config['prefix'], '') for x in fnames]
    middles = [x.replace('.fepout', '') for x in middles]
    # Remove first group of digits, then take all trailing lowercase letters
    nums, suffixes = set(), set()
    counter = 0 # Start at 'a'
    for s in middles:
        m = re.match(r'([0-9]+)([a-z]+)?', s)
        if m is None:
            print(f"generate_next_prefix: Filenames {fnames} don't follow the format I was expecting.")
            exit(1)
        matches = m.groups()
        nums.add(matches[0])
        if matches[1] is not None:
            suffixes.add(matches[1])

    # Start counting at the last caller-provided suffix
    suffixes = list(sorted(suffixes))
    if len(suffixes) > 0:
        counter = from_base26(suffixes[-1])
    highest_num_str = sorted(nums)[-1]

    # Pick a name without colliding with existing simulation output
    # Ensure we do at least one iteration
    out_prefix = fnames[0].replace('.fepout', '')
    # If there is already simulation output in the proposed filenames, then try another filename
    while any_isfile([f'{out_prefix}.fepout', f'{out_prefix}.restart.coor', f'{out_prefix}.restart.vel']):
        counter += 1
        out_prefix = f"{config['prefix']}{highest_num_str}{to_base26(counter)}"
    return out_prefix
        

def get_latest_restart_prefix(prefix):
    files = list(sorted(glob(f"{prefix}*.coor")))
    if len(files) == 0:
        return None
    return files[-1].replace('.coor', '')

# Returns the prefix of the restart files we should be using to start or continue this window.
# In order of highest to lowest priority:
#   1. The newest restart file from this window
#   2. The newest restart file from the previous window (if this window does not start a group)
#       - or -
#      The start_config as specified in the config YAML
def get_correct_restart_prefix(config, lambda_index):
    # This isn't the first lambda window, but it could still be starting a group
    # First lambda is by definition first window in group
    # Get newest restart file from this window
    this_prefix = get_latest_restart_prefix(f"{config['prefix']}{lambda_index:03d}")
    start_config_prefix = get_latest_restart_prefix(f"{config['start_config_prefix']}")
    last_window_prefix = get_latest_restart_prefix(f"{config['prefix']}{lambda_index - 1:03d}") if lambda_index > 0 else None
    is_first_window_in_group = lambda_index % config['windows_per_group'] == 0
    if this_prefix is not None:
        return this_prefix
    
    # No restart files for this window, so 
    if is_first_window_in_group:
        return start_config_prefix
    else:
        if last_window_prefix is None:
            print(f"Error: No restart files for previous lambda window index {lambda_index - 1}. That means I can't start this window.",
                file=sys.stderr)
            exit(1)
        return last_window_prefix      


def main():
    ap = argparse.ArgumentParser(description="""Dynamically generate a NAMD config file to start or continue an IDWS FEP window.

Sample configuration, in YAML format, that you would provide:

lambda_schedule: [0, 0.1, 0.2, 0.22, 0.25, 0.3, 0.35, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0]
prefix: myfep
start_config_prefix: prod1.restart
restart_frequency: 10000
total_steps_per_window: 2000000
windows_per_group: 5
per_group_equil_steps: 500000
per_window_equil_steps: 200000

We assume that you have a config_myfep.namd file that will be included in the generated NAMD configuration.""")
    ap.add_argument('config_yaml', help='YAML file describing this FEP calculation')
    # ap.add_argument('fepout_files', nargs='?', help='NAMD fepout files for a single window')
    ap.add_argument('lambda_index', type=int, help='Index in lambda_schedule for this window. We start counting at 0')
    args = ap.parse_args()

    # Read the config .yaml file
    # Look at the .fepout file we've got so far and decide how much more calculation we need to do
    with open(args.config_yaml) as f:
        config = safe_load(f)
    # print(config, file=sys.stderr)
    
    # Do a bunch of sanity checks as it is easy for the user to screw up
    if args.lambda_index < 0 or args.lambda_index >= len(config['lambda_schedule']):
        print(f"Error: Lambda index {args.lambda_index} is out of range - {args.config_yaml} only has {len(config['lambda_schedule'])} lambdas",
              file=sys.stderr)
        exit(1)

    if config['lambda_schedule'][0] != 0.0:
        print(f"Error: First lambda should be 0 and not {config['lambda_schedule'][0]}. I am finicky like that.", file=sys.stderr)
        exit(1)
        
    if config['lambda_schedule'][-1] != 1.0:
        print(f"Error: Last lambda should be 1.0 and not {config['lambda_schedule'][-1]}. Sorry.", file=sys.stderr)
        exit(1)
        
    if len(set(config['lambda_schedule'])) != len(config['lambda_schedule']):
        print("Error: Looks like there are duplicated lambda values. That's too weird for me.", file=sys.stderr)
        exit(1)

    # Determine what our lambdas are based on the provided index
    # We have already ensured that the index is something within range
    lambda1 = config['lambda_schedule'][args.lambda_index]
    
    # If we are already at the last lambda, which should be 1.0, then this is the backward IDWS last window
    # and we won't actually enable IDWS for this window
    if lambda1 == 1.0:
        lambda2 = config['lambda_schedule'][-2]
    else:
        lambda2 = config['lambda_schedule'][args.lambda_index+1]
   
    # No IDWS for the first window since there isn't any previous window
    # And none for the very last window since that's the backwards last window
    if args.lambda_index > 0 and lambda1 != 1.0:
        lambda_idws = config['lambda_schedule'][args.lambda_index-1]
    else:
        lambda_idws = None
        
    print(f'Lambdas are {lambda1} {lambda2} {lambda_idws}', file=sys.stderr)
        
    # How many steps in total have already been done for this window?
    total_steps_here = 0
    fepout_files = glob(f"{config['prefix']}{args.lambda_index:03d}*.fepout")
    for fepout_fname in fepout_files:
        steps_here = num_steps_in_fepout(fepout_fname)
        print(f'Num steps spanned by {fepout_fname}: {steps_here}', file=sys.stderr)
        total_steps_here += steps_here
        if steps_here < config['restart_frequency']:
            print(f"Warning: {fepout_fname} contains only {steps_here} steps, which is less than the restart frequency of {config['restart_frequency']}.",
                file=sys.stderr)

        this_lambda1, this_lambda2, this_lambda_idws = lambdas_for_fepout(fepout_fname)
        # We might not have any lambdas reported yet in the fepout if its job got killed really early,
        # so just take whatever lambdas specified in the config.yaml
        if this_lambda1 is None and steps_here > 0:
            this_lambda1, this_lambda2, this_lambda_idws = lambda1, lambda2, lambda_idws

        # Ensure that the lambdas are the same across all fepouts
        # If there is no lambda1 in this_fepout, assume that we haven't even reached the end of alchEquilSteps and that
        # the lambdas specified in the config.yaml are correct.
        if (lambda1, lambda2, lambda_idws) != (this_lambda1, this_lambda2, this_lambda_idws):
            print('Error: Different lambdas encountered within the same window among the config and existing fepout files.',
                file=sys.stderr)
            exit(1)
        lambda1, lambda2, lambda_idws = this_lambda1, this_lambda2, this_lambda_idws
   
    # Decide what we should continue from. Either the original configuration as specified
    # in the config YAML or the previous window
    input_prefix = get_correct_restart_prefix(config, args.lambda_index)
    if input_prefix is None:
        print("Error: Could not figure out coordinates to start from.", file=sys.stderr)
        print(f"Error: Please check that your {args.config_yaml} specifies files that actually exist.",
            file=sys.stderr)
        exit(1)

    total_steps_here -= total_steps_here % config['restart_frequency']
    total_steps_needed = config['total_steps_per_window'] - total_steps_here
    # Did we even make it past the alchEquilSteps? If not, do some more
    is_first_window_in_group = args.lambda_index % config['windows_per_group'] == 0
    equil_steps = config['per_group_equil_steps'] if is_first_window_in_group else config['per_window_equil_steps']
    # If we did some FEP steps but there's no restart file, we need to start from the beginning
    # Don't count the steps done since the last restart file
    equil_steps -= total_steps_here
    # Do any remaining equilibration steps
    equil_steps = 0 if equil_steps < 0 else equil_steps
    equil_steps_str = f'alchEquilSteps {equil_steps}' if equil_steps > 0 else ''
    # Add a minimization step at the start of each group if user asked for it
    minimize_str = f"minimize {config['minimize']}" if is_first_window_in_group and config['minimize'] > 0 else ''
    
    # If we need to 
    if total_steps_needed > 0:
        print(f'We need {total_steps_needed} more steps.', file=sys.stderr)
        # Ensure that the lambda1 specified in any existing fepouts is correct
        if config['lambda_schedule'].index(lambda1) != args.lambda_index:
            print(f'Error: lambda {lambda1} not in the right place in the configuration lambda_schedule', file=sys.stderr)
        common_namd = f"common_{config['prefix']}.namd"
        print(f'Including {common_namd} in this new NAMD config.', file=sys.stderr)
        print(f'Resuming from: {input_prefix}', file=sys.stderr)
        next_prefix = f"{config['prefix']}{args.lambda_index:03d}"
        # If any existing fepout files, don't overwrite them
        if len(fepout_files) > 0:
            next_prefix = generate_next_prefix(config, fepout_files)
        print('Next prefix is:', next_prefix, file=sys.stderr)
        # Generate a new NAMD config file from a supplied template
        # When generating a NAMD config where none existed before, keep track of window groups...usually every 5 windows
        idws_string = f"""alchLambdaIDWS {lambda_idws}""" if lambda_idws is not None else ''
        out_namd = f"""# FEP from {lambda1} to {lambda2} (IDWS {lambda_idws})
set inputname {input_prefix}
set outputname {next_prefix}
alchLambda {lambda1}
alchLambda2 {lambda2}
{idws_string}
source {common_namd}
# If the fepout frequency is higher than the restart frequency there could be duplicated steps
# if resuming.  You might consider removing any such duplicated steps in postprocessing.
firsttimestep {total_steps_here}
{minimize_str}
{equil_steps_str}
run {total_steps_needed}
"""
        output_namd_filename = f'{next_prefix}.namd'
        with open(output_namd_filename, 'w') as f:
            print(out_namd, file=f)
        print(f'Wrote new NAMD config to {output_namd_filename}.', file=sys.stderr)
        # Return the NAMD config filename we wrote on stdout so it can be used in a batch script
        print(output_namd_filename)

if __name__ == '__main__':
    main()
