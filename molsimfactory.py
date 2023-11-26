#!/usr/bin/env python3
#
# MD simulation parameter scanner, biased toward NAMD
#
# If you want to run the same simulation many times, varying one or more parameters each time,
# use this to generate the starting points in an organized fashion, each in an individual
# directory. This is done by providing a template starting point and specifying parameters
# to vary. 
#
# Here is an example of the YAML input:
#    title: SK2 N2O in selectivity filter ABF
#    parameters: {'num': {'start': 1, 'stop': 10, 'num': 10}}
#    psf: sm01_1n2o.psf
#    pdb: sm01_1n2o.pdb
#    restart_prefix: prod1n2opr2.restart
#    static_files: ['toppar']
#    template_files: ['pr3abf.colvars', 'pr3abf.namd']
#    submit_script: submit-abf-bender.sh
#    submit_command: sbatch
#    autoexec_commands: []
#
# Here is what each of those lines does:
#    title: Does nothing except being inserted in bits of output files
#    parameters: Specifies a parameter space to explore. Each tuple of parameters
#       will be given a separate generated run.
#    psf, pdb: Your actual molecular system for NAMD
#    restart_prefix: Starting point for your simulation. Should have coor, vel, xsc.
#    static_files: List of files and directories to symlink into each run directory
#    template_files: These are run through Jinja2 templating and copied into each
#       run directory. All the keys and values in this YAML file, as well as the
#       particular tuple of parameters, are available to the templating engine.
#    submit_script: A submit script you wrote for your particular batch queue
#       environment that will submit a single run. It will also be run through the
#       template engine and copied into its run directory.
#    submit_command: The command you use to submit something to the batch queue
#       system. This will be used in the generation of a "do-submit-all.sh" shell
#       script to submit all the runs.
#    autoexec_commands: List of commands to run for each generated run, during the
#       generation process, and not when the jobs are actually submitted.
#
# Tom Joseph, U. Penn
import argparse
import yaml
import itertools
import numpy as np
import os
import sys
import shutil
from jinja2 import Template, Environment, StrictUndefined, FileSystemLoader


def main():
    ap = argparse.ArgumentParser(description='Try a bunch of parameters for NAMD simulations')
    ap.add_argument('conf', help='Configuration file in YAML format describing how to generate runs')
    ap.add_argument('outdir', help='Name of directory to generate all these runs into')
    args = ap.parse_args()

    f = open(args.conf)
    conf = yaml.load(f, Loader=yaml.FullLoader)
    f.close()

    print(f"Title: {conf['title']}", file=sys.stderr)

    # Generate parameters for the runs, according to the user-supplied ranges
    space = generate_parameter_space(conf['parameters'])
    # Copy over relevant bits of the configuration file so they are available in the templates
    # Complete our list of static files that we will symlink
    if 'psf' in conf:
        space['psf'] = conf['psf']
        conf['static_files'].append(conf['psf'])
    if 'pdb' in conf:
        space['pdb'] = conf['pdb']
        conf['static_files'].append(conf['pdb'])

    space['title'] = conf['title']
    
    # Treat restart files separately, because user might choose not to use those
    if 'restart_prefix' in conf:
        space['restart_prefix'] = conf['restart_prefix']
        conf['static_files'].append(f"{conf['restart_prefix']}.coor")
        conf['static_files'].append(f"{conf['restart_prefix']}.vel")
        conf['static_files'].append(f"{conf['restart_prefix']}.xsc")

    # Keep a copy of the YAML file used to generate these runs
    conf['static_files'].append(args.conf)

    # Sanity check to see if all the files exist
    ensure_files_exist([*conf['static_files'], *conf['template_files']])

    # Do the cartesian product of lists of parameter values
    runs = expand_parameter_space(space)

    print(f"Generating {len(runs)} runs in {args.outdir}.", file=sys.stderr)
    # XXX: Sorry for the magic number. It was chosen because the format of run directories is
    # run%05d, and the idea was for the directory names to be easily sorted
    if len(runs) > 99999:
        print('Warning: You have a lot of runs, which means a ton of directories will be generated.',
            file=sys.stderr)
    os.makedirs(args.outdir, exist_ok=True)

    # runs is a list of dicts, one per run, that contain the parameter values
    conf['template_files'].append(conf['submit_script'])
    templates = load_templates(conf['template_files'])

    # Copy the static files into the top-level run directory.
    # This is done so that only that directory needs to be copied over to e.g. a cluster.
    gather_and_copy_static_files(conf['static_files'], args.outdir)

    # Now that we've read all the templates into memory, we can change the working
    # directory to our output directory, but saving the original path so we can
    # return to it
    os.chdir(args.outdir)
    out_wd = os.getcwd()

    # Get absolute paths of the static files we are going to symlink into each run directory.
    # We will symlink specifically to the copies we made, and not the original files.
    # This way, you only need to copy the top-level run directory specified by args.outdir,
    # and not its parent directory.
    static_files = [os.path.abspath(os.path.basename(fname)) for fname in conf['static_files']]

    # Generate the preamble to a master submission script
    # This automates the process of submitting all the runs to your favorite batch queuing system.
    master_submit_fname = 'do-submit-all.sh'
    master_submit_script = open(master_submit_fname, 'w')
    print('#!/bin/bash', file=master_submit_script)
    print(f"# {conf['title']}", file=master_submit_script)
    print('# Parameters:', file=master_submit_script)
    print(f"# {space}", file=master_submit_script)
    print(f"SUBMIT_COMMAND={conf['submit_command']}", file=master_submit_script)
    print(f"SUBMIT_SCRIPT={conf['submit_script']}", file=master_submit_script)

    # We use a counter to name the individual run directories
    run_counter = 1
    for run in runs:
        # Make the individual run directory, without caring if it already exists
        os.chdir(out_wd)
        run_wd = f"run{run_counter:05d}"
        os.makedirs(run_wd, exist_ok=True)
        os.chdir(run_wd)

        # Render templates like the NAMD config file
        run['run_counter'] = run_counter
        # Allow templates to access configuration parameters
        for k, v in conf.items():
            run[k] = v
            
        for fname, template in templates.items():
            with open(fname, 'w') as f:
                print(template.render(run), file=f)

        # Make symlinks to static files such as PSF, PDB, restart files.
        # We don't copy these because it would take up a lot of disk space.
        make_symlinks(static_files)

        # Make a YAML file that says what the parameters are
        with open('params.yaml', 'w') as f:
            print(yaml.dump(run), file=f)

        # Run user-specified commands in each run directory
        if 'autoexec_commands' in conf:
            for command in conf['autoexec_commands']:
                os.system(command)

        # Generate a line in the master submit script for this run
        print(f"cd {run_wd} ; $SUBMIT_COMMAND < $SUBMIT_SCRIPT ; cd ..", file=master_submit_script)
        run_counter += 1

    master_submit_script.close()
    print(f'Script to do all the job submissions: {args.outdir}/{master_submit_fname}')

# Iterate through dict and expand each term as necessary.
# For example:
#   {'a': {'start': 1, 'stop': 3, 'num': 3}, 'b': ['foo', 'bar']}
#       ==>
#   {'a': [1, 2, 3], 'b', ['foo', 'bar']}
# This is suitable for passing into expand_parameter_space which will make the
# cartesian product representing the individual runs.
def generate_parameter_space(parameter_description):
    out = {}
    for k, v in parameter_description.items():
        if type(v) is dict:
            # Expand out a range (really a linear space) of values for this parameter, if
            # that's what the user asked for by specifying start, stop, and num parameters
            if set(['start', 'stop', 'num']).issubset(set(v)):
                out[k] = np.linspace(v['start'], v['stop'], v['num'], endpoint=True).tolist()
            else: # We are duck-typing these
                out[k] = v
        elif type(v) in (float, int, str, list):
            out[k] = v
        else:
            print(f'I have no idea what is going on with these parameters: {k}: {v}')

    return out


# Make cartesian product of a parameter space, suitable for generating runs.
# For example:
#   {'a': [1, 2, 3], 'b', ['foo', 'bar']}
#       ==>
#   [{'a': 1, 'b': 'foo'},
#    {'a': 2, 'b': 'foo'},
#    {'a': 3, 'b': 'foo'},
#    {'a': 1, 'b': 'bar'},
#    {'a': 2, 'b': 'bar'},
#    {'a': 3, 'b': 'bar'}]
def expand_parameter_space(space):
    # Unzip the parameter space dict
    #   {'foo': [1,2,3], 'bar':['A','B','C']} -> ['foo', 'bar'], [1,2,3], ['A','B',C']
    # Ordering must be preserved or we won't be able to match up parameters to their values later
    names, values = [], []
    for k, v in space.items():
        names.append(k)
        if type(v) is not list:
            values.append([v,])
        else:
            values.append(v)

    # Calculate the cartesian product of all those value lists
    expanded = itertools.product(*values)

    # Annotate each entry of the cartesian product with the parameter names so that we can
    # pass them to the templating engine
    runs = []
    for run_values in expanded:
        run = {}
        for i in range(len(names)):
            run[names[i]] = run_values[i]
        runs.append(run)
    return runs

# Convenience function to check if a list of files exist, and warn if any don't  exist
def ensure_files_exist(fnames):
    for fname in fnames:
        if os.path.exists(fname) is False:
            print(f'Warning: {fname} does not seem to exist, which may be problematic.', file=sys.stderr)

# Load template files and instantiate Jinja2 templates.
# This is a convenience function
def load_templates(fnames):
    env = Environment(loader=FileSystemLoader('.'), undefined=StrictUndefined)
    templates = {}
    for fname in fnames:
        templates[fname] = env.get_template(fname)
    return templates


# Make symlinks of files into the current directory.
# If a symlink already exists, such as when we are regenerating
# into the same directory, remove it in case it points to the wrong place
def make_symlinks(fnames):
    for fname in fnames:
        try:
            os.remove(os.path.basename(fname))
        except:
            pass
        os.symlink(os.path.relpath(fname), os.path.basename(fname))


# Gathers and copies static files and directories into one place, overwriting what's there.
# It is this complicated because it must handle copying both directory trees and individual
# files, and dealing with errors when the destination path exists.
def gather_and_copy_static_files(fnames, outdir):
    for fname in fnames:
        dest = os.path.join(outdir, os.path.basename(fname))
        if os.path.exists(dest):
            try:
                shutil.rmtree(dest)
            except NotADirectoryError:
                os.remove(dest)
        try:
            shutil.copytree(fname, dest)
        except NotADirectoryError:
            shutil.copy(fname, dest)


if __name__ == '__main__':
    main()
