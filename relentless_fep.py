#!/usr/bin/env python3
"""Relentless FEP

Eases administrative burden of large FEP calculations by automating, to some degree, the
resumption of interrupted NAMD runs.

Given a description of the overall FEP calculation and a particular window's .fepout files,
generates a new .namd file that starts or continues the calculation for this window.
Intended to tolerate unexpectedly interrupted calculations, such as when your jobs are
preempted.

"Quick" start
-------------

This script is only one part of a FEP calculation infrastructure. All it does is examine
the output of a single FEP window run (if it exists), then generates a NAMD configuration
file that, when run to completion, will finish that window. You are expected to run the
windows in order yourself, which usually means writing a batch job submission script.

The FEP calculation is divided into groups of windows. Each group begins from the same
initial simulation state. This allows for better parallelization at the cost of more
overall computation required due to an increased need for equilibration. For example,
if your FEP calculation needs 120 windows, running all the windows in sequence will take
a long time. But if you can run 0-4, 5-9, 10-14, etc. simultaneously on your cluster, you
can save a lot of wallclock time.

Of course, if you want to run all the windows sequentially, you can just say that each
group should have 121 or more windows. (That last window is the final IDWS window in the
reverse direction.)

Each FEP calculation requires:

  1. A configuration file in YAML format listing your requirements for the entire
     calculation
  2. A NAMD configuration that will be included in each generated NAMD configuration file
  3. Associated files such as PSF/PDB, parameter files, colvars files, and whatever else
     your particular simulations need

You must choose a prefix string that identifies your YAML and NAMD config files. For
example, if you choose the prefix 'myfep', your YAML file must be named 'myfep.yaml' and 
your NAMD config 'myfep.include.namd'. The output files will have the form 'myfep000.fepout',
'myfep000.log', 'myfep000.colvars.traj', etc.

Here is an example YAML file, with 21 total windows (including the final reverse windows),
split into groups of 5:

```
lambda_schedule: [0, 0.1, 0.12, 0.17, 0.2, 0.22, 0.25, 0.3, 0.35, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 1.0]
start_config_prefix: mineq
restart_frequency: 10000
total_steps_per_window: 1500000
windows_per_group: 5
minimize: 100
per_group_equil_steps: 500000
per_window_equil_steps: 200000
```

And a description of the parameters:
  - "lambda_schedule" must span 0 to 1.0.
  - "start_config_prefix" identifies the simulation input to begin from with the first window in
    each group. Here this means mineq.coor, mineq.vel, mineq.xsc. You will specify the PSF file
    in the NAMD config myfep.include.namd.
  - "restart_frequency" is the NAMD restart frequency that you specified in the NAMD config. This
    is important so relentless_fep can determine exactly how many steps are necessary to run
    when a simulation is interrupted and resumed.
  - "total_steps_per_window" is the goal number of simulation steps per FEP window. In this case,
    assuming a 2 fs timestep, each window would be 3 ns. This includes both equilibration and
    production simulation.
  - "windows_per_group" specifies how many windows each group should contain. A group of windows
    can be run in parallel with the other groups of windows, because each group starts from
    the same simulation configuration (here, mineq.coor, mineq.vel, mineq.xsc).
  - "minimize" allows for minimization of the simulation prior to dynamics for only the first
    window in each group. This is useful if you edit the simulation input system by hand, such as
    to modify a functional group.
  - "per_group_equil_steps" is the number of equilibration steps to run at the start of the
    group. This is *not* in addition to the per_window_equil_steps. This is the value passed to
    NAMD using the alchEquilSteps option.
  - "per_window_equil_steps" is the number of equilibration steps for all windows except the first
    one in a group.

When relentless_fep.py is run for window 0 as follows:

```
$ python3 relentless_fep.py myfep 0
```

...it will generate a NAMD config file named myfep000.namd that looks like this:

```
# FEP from 0.0 to 0.1 (IDWS None)
set inputname mineq
set outputname myfep000
alchLambda 0.0
alchLambda2 0.1
source myfep.include.namd
firsttimestep 0
minimize 100
alchEquilSteps 500000
run 1500000
```

Your myfep.include.namd should include everything else NAMD needs to run your simulation.

Let's say you ran this simulation, but then you had a power outage after 400000 steps. All you
now must do is re-run relentless_fep.py in exactly the same way as before:

```
$ python3 relentless_fep.py myfep 0
```

It will generate a NAMD config file named myfep000a.namd that looks like this:

```
# FEP from 0.0 to 0.1 (IDWS None)
set inputname myfep000.restart
set outputname myfep000a
alchLambda 0.0
alchLambda2 0.1
source myfep.include.namd
firsttimestep 400000
alchEquilSteps 100000
run 1100000
```

Naturally, it is cumbersome to run these simulations by hand. So, you should write a batch
submission script that automates the process. Since compute clusters vary in how you must
interact with them, this is unfortunately highly dependent on your local setup. But, as an
example, here is the skeleton of such a script that would work with the SLURM batch queue:

```
#!/bin/bash
NAMD=$HOME/bin/namd2
PREFIX=myfep
FIRST=0
LAST=20
GROUP_SIZE=5
for group_start in `seq $FIRST $GROUP_SIZE $LAST`; do
    group_end=$(($group_start + $GROUP_SIZE - 1))
    if (( $group_end > $LAST )); then
        $group_end = $LAST
    fi

    SCRIPT() { cat <<EOF
#!/bin/bash
#SBATCH -J ${PREFIX}_${group_start}_${group_end}
for m in \`seq $group_start $group_end\`; do
    m=\`printf %03d \$m`
    namdconf=\$(python3 relentless_fep.py ${PREFIX}.yaml \${m})
    # If relentless_fep did not generate a NAMD config file, it's done with this window,
    # so don't try to run it
    if [ -z \$namdconf ]; then continue; fi
    namdlog=\$(basename \$namdconf .namd).log
    $NAMD \${namdconf} > \${namdlog}
done
EOF
    }

    # Submit the job
    SCRIPT | sbatch
done
```

If a calcuation fails, you can just resubmit the script without any other changes. Often
the batch queue system can be configured to do this automatically for you.


"""
import argparse
import sys
import re
from yaml import safe_load
import os.path
from glob import glob
from string import ascii_lowercase
from functools import reduce

# How many MD steps are represented by a given .fepout file?
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

    # NAMD will generate a fepout file with no steps in it, apparently?
    if len(a_few_steps) < 2:
        total_steps = 0
    else:
        # Avoid fencepost error because the FepEnergy for the first timestep is not recorded
        total_steps = max_step_num - min_step_num + (a_few_steps[1] - a_few_steps[0])
    return total_steps


# Identifies the lambda values represented in a .fepout file
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


# Returns True if specified fepout file contains "#Free energy change for lambda window" line.
# We use this to ensure that the fepout file set for a given lambda value is actually considered
# complete by a fepout parser.
def fepout_contains_footer(fname):
    with open(fname) as f:
        for line in f:
            if line.startswith('#Free energy change for lambda window'):
                return True
    return False


# Converts a number to a series of letters in base 26, like MS Excel does to label columns
# These are used for suffixes to generated namd config files rather than the numbers
# themselves, because that would be confusing.
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
        

# Returns the newest restart file of a given prefix.
# Assumes that you are using base26 suffixes, so that an alphabetical sort works.
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

Sample configuration myfep.yaml, in YAML format, that you would provide:

lambda_schedule: [0, 0.1, 0.2, 0.22, 0.25, 0.3, 0.35, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0]
start_config_prefix: prod1
restart_frequency: 10000
total_steps_per_window: 2000000
windows_per_group: 5
per_group_equil_steps: 500000
per_window_equil_steps: 200000

You must have a myfep.include.namd file that will be included in the generated NAMD configuration.""",
        formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument('prefix', help='Prefix identifying the files for this FEP calculation')
    # ap.add_argument('fepout_files', nargs='?', help='NAMD fepout files for a single window')
    ap.add_argument('lambda_index', type=int, help='Index in lambda_schedule for this window. We start counting at 0')
    args = ap.parse_args()

    # Read the config .yaml file
    # Look at the .fepout file we've got so far and decide how much more calculation we need to do
    config_yaml = f'{args.prefix}.yaml'
    config['prefix'] = args.prefix
    common_namd = f"{config['prefix']}.include.namd"

    # Ensure the config YAML and NAMD config files exist
    if os.path.exists(config_yaml) is False:
        print(f'Error: Configuration file {config_yaml} not found.', file=sys.stderr)
    if os.path.exists(common_namd) is False:
        print(f'Error: NAMD configuration file {common_namd} not found.', file=sys.stderr)

    with open(config_yaml) as f:
        config = safe_load(f)    

    # Do a bunch of sanity checks as it is easy for the user to screw up
    if args.lambda_index < 0 or args.lambda_index >= len(config['lambda_schedule']):
        print(f"Error: Lambda index {args.lambda_index} is out of range - {args.config_yaml} only has {len(config['lambda_schedule'])} lambdas",
              file=sys.stderr)
        exit(1)

    endpoint_lambdas = config['lambda_schedule'][0], config['lambda_schedule'][-1]
    if endpoint_lambdas != (0.0, 1.0) and endpoint_lambdas != (1.0, 0.0):
        print(f"Error: Lambda schedule should begin and end with 0.0 and 1.0, or vice versa.", file=sys.stderr)
        print(f"       They are actually {endpoint_lambdas}. I am finicky like that.", file=sys.stderr)
        exit(1)
       
    if len(set(config['lambda_schedule'])) != len(config['lambda_schedule']):
        print("Error: Looks like there are duplicated lambda values. That's too weird for me.", file=sys.stderr)
        exit(1)

    # Determine what our lambdas are based on the provided index
    # We have already ensured that the index is something within range
    lambda1 = config['lambda_schedule'][args.lambda_index]
    
    # If we are already at the last lambda, then this is the backward IDWS last window
    # and we won't actually enable IDWS for this window
    if args.lambda_index == len(config['lambda_schedule']) - 1:
        lambda2 = config['lambda_schedule'][-2]
    else:
        lambda2 = config['lambda_schedule'][args.lambda_index+1]
   
    # No IDWS for the first window since there isn't any previous window
    # And none for the very last window since that's the backwards last window
    lambda_idws = config['lambda_schedule'][args.lambda_index-1] if lambda1 not in (0.0, 1.0) else None
        
    print(f'Lambdas are {lambda1} {lambda2} {lambda_idws}', file=sys.stderr)
        
    # How many steps in total have already been done for this window?
    # Also, if we have run all the steps user asked for, check to see whether the
    # final fepout file contains the footer. If it doesn't, generate a config file
    # run NAMD one last time to generate the fepout.
    total_steps_here = 0
    footer_present = False
    fepout_files = glob(f"{config['prefix']}{args.lambda_index:03d}*.fepout")
    for fepout_fname in fepout_files:
        steps_here = num_steps_in_fepout(fepout_fname)
        if fepout_contains_footer(fepout_fname):
            footer_present = True

        print(f'Num steps spanned by {fepout_fname}: {steps_here}', file=sys.stderr)
        total_steps_here += steps_here
        if steps_here < config['restart_frequency']:
            print(f"Warning: {fepout_fname} contains only {steps_here} steps, which is less than the restart frequency of {config['restart_frequency']}.",
                file=sys.stderr)
        # Handle the case where the previous fepout file is empty entirely.
        # This can happen with AWS spot instances that are extremely short lived.
        # Fall through.
        if steps_here == 0:
            break

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
            print(f'Previous lambda1, lambda2, lambda_idws = {lambda1}, {lambda2}, {lambda_idws}', file=sys.stderr)
            print(f'This lambda1, lambda2, lambda_idws = {this_lambda1}, {this_lambda2}, {this_lambda_idws}', file=sys.stderr)
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
    # But don't minimize if we are resuming within the same window
    if input_prefix == config['prefix']:
        minimize_str = ''
    
    # Perhaps we have enough steps but still need to generate a footer
    if total_steps_needed <= 0 and footer_present is False:
        print('There is no footer line in this set of fepouts, so we will run slightly longer.', file=sys.stderr)
        total_steps_needed = 100

    # If we need to run more steps, do it
    if total_steps_needed > 0:
        print(f'We need {total_steps_needed} more steps.', file=sys.stderr)
        # Ensure that the lambda1 specified in any existing fepouts is correct
        if config['lambda_schedule'].index(lambda1) != args.lambda_index:
            print(f'Error: lambda {lambda1} not in the right place in the configuration lambda_schedule', file=sys.stderr)
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
