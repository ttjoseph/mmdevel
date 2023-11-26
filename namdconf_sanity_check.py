#!/usr/bin/env python
#
# NAMD configuration file sanity checks
#
# This is pretty much a transliteration of Ezry's quality control scripts
# from awk into Python, except we use the Tcl intepreter to parse the
# NAMD config file.
#
# Tom Joseph, U. Penn
import sys
import argparse
import tkinter
import re
import pathlib
import numpy as np
import warnings
import MDAnalysis as mda

# MDAnalysis will complain about large PDB files. It is correct to do so, but our
# workflow makes them unavoidable
warnings.filterwarnings("ignore", category=UserWarning, module='MDAnalysis')

def main():
    ap = argparse.ArgumentParser(description='Sanity checker for NAMD FEP MD files')
    ap.add_argument('filename', nargs='+')
    args = ap.parse_args()

    for fname in args.filename:
        print(f"--- Sanity check: {fname} ---")
        load_namdconf(fname)


# We automate the extraction of values from the NAMD config file because it is pretty much
# just a Tcl script. So we run it through the Tcl interpreter, defining all commands to just
# save the values in a Tcl array which we then extract and examine in the friendly confines
# of Python.
TCL_INIT_SCRIPT = """
rename unknown unknown_orig
proc unknown {args} {
    global OUT
    set cmd [lindex $args 0]
    set cmd_args [lrange $args 1 end]
    lappend OUT($cmd) $cmd_args
}
"""

def load_namdconf(fname):
    with open(fname, 'r') as f:
        namdconf_str = f.read()

    # Spawn a Tcl interpreter instance and run our script in it
    # to prepare it to ingest a NAMD config file
    tcl = tkinter.Tcl()
    tcl.eval(TCL_INIT_SCRIPT)
    
    # Here is where we actually run the NAMD config file.
    tcl.eval(namdconf_str)
    
    # Extract data from the Tcl environment in this absurd manner
    # It only gives us a string rather than a more structured format,
    # so we need to do a little more parsing
    output = tcl.eval('array names OUT')
    conf = {}
    for k in output.split():
        data = tcl.eval(f"array get OUT {k}")
        # Strip any enclosing braces if the value returned a Tcl list, which
        # is really a string with some syntactic sugar
        # TODO: So this is wrong if we are extracting a Tcl list
        # TODO: there has got to be a better way of getting data out of Tcl
        data = data[1:] if data.startswith('{') else data
        data = data[:-1] if data.endswith('}') else data
        tokens = re.sub('#.*', '', data).strip().split()
        # NAMD allows case insensitivity for parameters
        conf[tokens[0].lower()] = ' '.join(tokens[1:])

    # Report some parameters to the user
    # Note that all keys in the conf dict are lowercase, no matter what they were
    # in the NAMD config file we parsed!
    if 'colvars' in conf and is_tcltrue(conf['colvars']):
        print("INFO: Uses colvars")
        colvarsconfig = pathlib.Path(conf['colvarsconfig'])
        if colvarsconfig.exists() is False:
            print('WARNING: But the colvars config file does not appear to exist')

    print(f"INFO: Runs for {conf['run']} steps")
    if is_tcltrue(conf['langevin']):
        print(f"INFO: Uses Langevin dynamics")
        print(f"INFO: Langevin period = {conf['langevinpistonperiod']}")
        print(f"INFO: Langevin decay = {conf['langevinpistondecay']}")

    # Nobody likes wrapAll
    if is_tcltrue(conf['wrapall']):
        print('ERROR: wrapAll is turned on. This is always a bad idea.')

    # If we're doing alchemy, sanity-check that
    if 'alch' in conf and is_tcltrue(conf['alch']):
        alchfile = pathlib.Path(conf['alchfile'])
        if alchfile.exists():
            print(f'INFO: FEP alchFile {conf["alchfile"]} exists')
            # Now count the atoms undergoing perturbation
            fep_u = mda.Universe(alchfile, topology_format='PDB')
            num_pert_atoms = 0
            if conf['alchcol'] == 'B':
                num_pert_atoms = np.count_nonzero(fep_u.atoms.tempfactors)
            elif conf['alchcol'] == 'O':
                num_pert_atoms = np.count_nonzero(fep_u.atoms.occupancies)
            else:
                print(f'ERROR: I dunno what alchCol {conf["alchcol"]} means')
            print(f'INFO: FEP number of atoms perturbed: {num_pert_atoms}')
            if num_pert_atoms > 100:
                print('WARNING: Perturbing {num_per_atoms} atoms is risky business, that is a lot!')
        else:
            print(f'ERROR: FEP alchFile {conf["alchfile"]} does not exist')

        if 'alchoutfreq' in conf and 'fullelectfrequency' in conf and 'nonbondedfreq' in conf:
            fullelectfreq, nonbondedfreq, alchoutfreq = int(conf['fullelectfrequency']), int(conf['nonbondedfreq']), int(conf['alchoutfreq'])
            if int(alchoutfreq / fullelectfrequency) != (alchoutfreq / fullelectfrequency):
                print(f'WARNING: alchOutFreq is not a multiple of fullElectFrequency')
                print(f'         See: https://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2020-2021/1487.html')
            elif int(alchoutfreq / nonbondedfreq) != (alchoutfreq / nonbondedfreq):
                print(f'WARNING: alchOutFreq is not a multiple of nonbondedFreq')
                print(f'         See: https://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2020-2021/1487.html')



# Many ways to say true in NAMD config files
def is_tcltrue(val):
    return val.strip().lower() in ['on', 'true', 'yes']

if __name__ == '__main__':
    main()


