#!/usr/bin/env python
import sys
import argparse
import csv
import re
import yaml
import numpy as np
import matplotlib.pyplot as plt
import cycler
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from ffparam_charges_vs_qm import load_prm, COULOMB


def calc_mm_interaction_energy(u, ligand_spec, prm, solute_spec, vdw_shift=0.0, patches=None):
    """
    Given an MDAnalysis universe, calculate the pairwise interaction energies between the specified
    ligand and the solute.
    :param u: MDAnalysis universe
    :param ligand_spec: Specification to find the ligand molecule, such as "resname APM"
    :param prm: Dict containing Lennard-Jones parameters by atom type, as read from ffparam_charges_vs_qm.load_prm
    :param solute_spec: Specification to find the solute, such as "protein or nucleic"
    :param vdw_shift: In L-J calculations only, increase interparticle distance by this many Angstroms
        for the sake of softcore FEP (see the documentation for alchVdwShiftCoeff in NAMD)
    :param patches: A dict containing modifications to the ligand to be applied, such as nudging charges.
        Should be like {SKE_N1_chargedelta: 0.1}, which increases the charge of the N1 atom of SKE by 0.1.
    :return: A tuple containing electrostatic and Lennard-Jones energies, as NxM arrays, with N ligand atoms and M frames
    """
    ligand_atoms = u.select_atoms(ligand_spec)
    solute_atoms = u.select_atoms(solute_spec)
    if len(ligand_atoms) == 0:
        print >>sys.stderr, 'Cannot find any ligand atoms. Did you provide the correct residue name?'
        sys.exit(1)
    if len(solute_atoms) == 0:
        print >>sys.stderr, 'Cannot find any solute atoms. Did you provide the correct spec?'
        sys.exit(1)

    ligand_rmin2 = np.array([prm[a.type]['rmin2'] for a in ligand_atoms])
    ligand_epsilon = np.array([prm[a.type]['epsilon'] for a in ligand_atoms])
    solute_rmin2 = np.array([prm[a.type]['rmin2'] for a in solute_atoms])
    solute_epsilon = np.array([prm[a.type]['epsilon'] for a in solute_atoms])

    # We will be precomputing the pairwise charge products, which we will later scale by inverse distance,
    # as well as part of the Lennard-Jones calculation.
    # We can do this otherwise horrible O(n^2) calculation because we know we have a relatively
    # small number of ligand atoms, so the size of the resulting matrix is not cataclysmic.
    # There has got to be a Numpy-array-based way to do this, but since we only do it once
    # we can tolerate the slowness.
    charge_prod = np.zeros((len(ligand_atoms), len(solute_atoms)))
    lj6_partial, eps_partial = np.zeros_like(charge_prod), np.zeros_like(charge_prod)
    # Fill in default of no patches
    if patches is None: patches = {}

    # Modify solute atom parameters according to user-specified patches.
    # This mechanism is different from the ligand one because in the ligand case we want to
    # modify all ligand molecules in the same way. In this case, we really just want one
    # atom to change. Obviously this means there is no clear way to only modify one ligand
    # molecule, but I guess that can be added later if needed.
    for solute_i in range(len(solute_atoms)):
        solute_atom = solute_atoms[solute_i]
        # HSD_123_HB2_chargedelta: -0.5
        chargedelta_key = '%s_%d_%s_chargedelta' % (solute_atom.resname, solute_atom.resid, solute_atom.name)
        if chargedelta_key in patches:
            solute_atom.charge += patches[chargedelta_key]

    # Modify ligand atom charges according to the patches specified by the user
    for ligand_i in range(len(ligand_atoms)):
        ligand_atom = ligand_atoms[ligand_i]
        # Modify charge of this ligand atom if user asked us to
        chargedelta_key = '%s_%s_chargedelta' % (ligand_atom.resname, ligand_atom.name)
        if chargedelta_key in patches:
            ligand_atom.charge += patches[chargedelta_key]

        ligand_charge = COULOMB * ligand_atom.charge
        l_rmin2, l_eps = ligand_rmin2[ligand_i], ligand_epsilon[ligand_i]
        for solute_i in range(len(solute_atoms)):
            charge_prod[ligand_i, solute_i] = ligand_charge * solute_atoms[solute_i].charge
            lj6_partial[ligand_i, solute_i] = (l_rmin2 + solute_rmin2[solute_i])**6
            eps_partial[ligand_i, solute_i] = (l_eps * solute_epsilon[solute_i])**0.5

    electro, lj = [], []

    for ts in u.trajectory:
        # Tell the user every so often which frame we're on
        if (ts.frame > 0 and ts.frame % 25 == 24) or ts.frame+1 == len(u.trajectory):
            print >>sys.stderr, 'Processing frame %d of %d.' % (ts.frame+1, len(u.trajectory))
        # Calculate pairwise distances, which thankfully we don't need to do manually, because
        # that would be incredibly slow in pure Python.
        dists = distance_array(ligand_atoms.positions, solute_atoms.positions)
        # ...and scale the precomputed charge products by those distances, summing over all solute
        # atoms to yield a per-ligand-atom electrostatic interaction energy.
        frame_electro = np.sum(charge_prod / dists, 1)
        # Then do the analogous thing for Lennard-Jones energies.
        # Lennard-Jones energy calculation was:
        # lj6 = ((a_rmin2 + prm[b.type]['rmin2']) / dist)**6
        # ligand_lj[ligand_i] += (a_epsilon * prm[b.type]['epsilon'])**0.5 * (lj6**2 - 2*lj6)

        # For FEP softcore purposes, shift the interparticle distance
        if vdw_shift != 0.0:
            dists += vdw_shift
        lj6 = lj6_partial / (dists**6)
        frame_lj = np.sum(eps_partial * (lj6**2 - 2*lj6), 1)

        electro.append(frame_electro)
        lj.append(frame_lj)

    return electro, lj


def scale_factors_for_lambda(lambda_val, alchElecLambdaStart=0.5, alchVdwLambdaEnd=0.5, annihilate=True):
    """
    Returns scale factors for electrostatic and L-J energies for given softcore FEP parameters and lambda

    :param lambda_val: The FEP lambda, which should be from 0.0 to 1.0
    :param alchElecLambdaStart: Same as with NAMD
    :param alchVdwLambdaEnd: Same with NAMD
    :return: electro_scale, lj_scale
    """

    # TODO: Distinguish between annihilation and exnihilation, because the meanings of the parameters
    # TODO: depend on it.
    # TODO: When annihilating, alchElecLambdaStart actually acts like a "alchElecLambdaEnd"

    if annihilate is True:
        lambda_val = 1.0 - lambda_val

    electro_scale, lj_scale = 0.0, 1.0
    if lambda_val >= alchElecLambdaStart:
        electro_scale = (lambda_val - alchElecLambdaStart) / (1.0 - alchElecLambdaStart)

    if lambda_val <= alchVdwLambdaEnd:
        lj_scale = lambda_val / alchVdwLambdaEnd

    return electro_scale, lj_scale


# TODO: Since this is just scaling, we don't have to sum across ligand atoms...do we?
def do_fep_scaling(electro, lj, lambda0, lambda1, alchElecLambdaStart, alchVdwLambdaEnd):
    electro_scale0, lj_scale0 = scale_factors_for_lambda(args.lambda0, alchElecLambdaStart, alchVdwLambdaEnd)
    electro_scale1, lj_scale1 = scale_factors_for_lambda(args.lambda1, alchElecLambdaStart, alchVdwLambdaEnd)

    # Now, scale the energies
    electro0, lj0, electro1, lj1 = [], [], [], []
    for i in range(len(electro)):
        electro0.append(electro[i] * electro_scale0)
        lj0.append(lj[i] * lj_scale0)
        electro1.append(electro[i] * electro_scale1)
        lj1.append(lj[i] * lj_scale1)

    return electro0, lj0, electro1, lj1


# TODO: If weights are extremely nonuniform, complain that the energies1 do not cover a similar enough
# TODO: part of the phase space represented by energies0, and therefore the reweighted energy is probably
# TODO: far off.
def reweight_energies(energies0, energies1):
    e0 = np.array(energies0)
    e1 = np.array(energies1)
    assert len(e0) == len(e1)

    beta = 1.0/0.6
    weights = np.exp(-beta * (e1 - e0))
    weights /= np.sum(weights)

    print >>sys.stderr, 'Then: %.3f  Now: %.3f' % (np.sum(e0)/np.size(e0), np.sum(e1 * weights))

    return(e1 * weights)


# TODO: This should just take a CSV filename or something
def generate_plot(datatype, data, dont_plot_atomtypes=''):
    # Make a list of labels specified by the user not to include in the graph, because they are boring
    # atoms like aliphatic hydrogens with a charge of +0.09
    dont_plot_labels = []
    for l in labels:
        for atomtype in dont_plot_atomtypes.split(','):
            if l.startswith('%s ' % atomtype):
                dont_plot_labels.append(l)
    print >>sys.stderr, 'Will not plot:', ', '.join(dont_plot_labels)

    # Try not to repeat line styles
    pretty_cycler = cycler.cycler(lw=[0.3, 0.8, 1.4]) * cycler.cycler('ls', ['-', ':']) * plt.rcParams['axes.prop_cycle']
    plt.rc('axes', prop_cycle=pretty_cycler)

    orig_d = pd.read_csv('%s_%s.csv' % (datatype, ligand_resname))
    print >> sys.stderr, 'Average %s energy: %.2f kcal/mol' % (datatype, np.sum(np.mean(orig_d)))

    # Generate plot
    d = orig_d.drop(dont_plot_labels, axis=1)
    plt.figure(figsize=(8, 4.5))  # 16:9 aspect ratio
    plt.autoscale(tight=True)
    plt.plot(d)
    plt.figtext(0.01, 0.99, 'Average %s energy, including any atom types not plotted: %.2f kcal/mol' % \
                (datatype, np.sum(np.mean(orig_d))), fontsize=5, verticalalignment='top')
    plt.title('%s energy for %s with %s (%s)' % (datatype.capitalize(), args.ligand_resname, args.solute_spec,
                                                 args.namdconf))
    plt.xlabel('Trajectory frame')
    plt.ylabel('Energy (kcal/mol)')

    plt.legend(list(d), fontsize=4, loc='lower left', bbox_to_anchor=(1.0, 0.0))
    plt.savefig('%s.pdf' % prefix)
    print >> sys.stderr, 'Wrote a pretty picture to %s.pdf' % prefix


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Calculates single point energies with MM, with respect to a ligand.')
    ap.add_argument('ligand_resname', help='Ligand residue name (e.g. SKE, APM)')
    ap.add_argument('namdconf', help='NAMD configuration containing "parameters" and "structure" keywords from which to glean PSF and force field filenames')
    ap.add_argument('dcd', nargs='+', help='DCD trajectories to go with the PSF file specified in namdconf')
    ap.add_argument('--solute-spec', default='protein or nucleic', help='Specify which part of the system we care about ligand interactions with')
    ap.add_argument('--dont-plot-atomtypes', default='HA1,HA3,HP,HGA2,HGAAM1', help='Atom types not to plot, separated by commas')
    ap.add_argument('--patches', help='YAML file containing patches to ligand parameters, such as charges')
    ap.add_argument('--lambda0', type=float, help='FEP lambda0 for the provided trajectory')
    ap.add_argument('--lambda1', type=float, help='FEP lambda1 for the provided trajectory')
    ap.add_argument('--alchElecLambdaStart', type=float, default=0.5, help='Governs lambda scaling, as it does in NAMD')
    ap.add_argument('--alchVdwLambdaEnd', type=float, default=0.5, help='Governs lambda scaling, as it does in NAMD')
    args = ap.parse_args()

    psf_filename, prm_filenames = None, []

    # Glean PSF and force field parameters from NAMD simulation input file
    print >>sys.stderr, 'Reading NAMD simulation config', args.namdconf
    with open(args.namdconf) as f:
        psf_re = re.compile(r'structure\s+([\w\-/.]+)')
        prm_re = re.compile(r'parameters\s+([\w\-/.]+)')
        for line in f:
            m = psf_re.match(line.strip())
            if m: psf_filename = m.group(1)
            m = prm_re.match(line.strip())
            if m: prm_filenames.append(m.group(1))

    print >>sys.stderr, 'Guessed PSF filename:', psf_filename
    print >>sys.stderr, 'Guessed force field parameter files:', ', '.join(prm_filenames)
    print >>sys.stderr, 'Solute spec:', args.solute_spec

    u = mda.Universe(psf_filename, args.dcd)

    prm = {}
    for prm_fname in prm_filenames:
        prm.update(load_prm(prm_fname))

    print >>sys.stderr, 'Loaded %d MM parameters.' % len(prm)

    patches = None
    if args.patches is not None:
        patches = yaml.load(open(args.patches))
        print >>sys.stderr, 'Patches:', patches

    ligand_spec = 'resname %s' % args.ligand_resname

    # Calculate the energies assuming no FEP and no patches
    print >>sys.stderr, 'Calculating energies with no patches.'
    electro_orig, lj_orig = calc_mm_interaction_energy(u, ligand_spec, prm, solute_spec=args.solute_spec)
    all_data = [('electro', electro_orig), ('lj', lj_orig)]

    # If we have patches, recalculate the energies
    if patches is not None:
        print >> sys.stderr, 'Calculating energies WITH patches.'
        electro_patched, lj_patched = calc_mm_interaction_energy(u, ligand_spec, prm, solute_spec=args.solute_spec,
                                                                 patches=patches)
        # Now we can reweight the patched energies
        electro_reweighted = reweight_energies([np.sum(en) for en in electro_orig],
                                               [np.sum(en) for en in electro_patched])
        # We actually don't need to do this, yet, because we don't change any L-J parameters
        lj_reweighted = reweight_energies([np.sum(en) for en in lj_orig],
                                          [np.sum(en) for en in lj_patched])

    # We just calculated the unscaled electrostatic and Lennard-Jones energies.
    # If this was a FEP run, we now need to scale the energies according to the lambda value.
    # Because electro and lj contain energies assuming there is no scaling (e.g. lambda is 0 for decoupling
    # or 1 for coupling)

    # DEBUG: show how the lambdas would scale
    # for l in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
    #     print >>sys.stderr, l, scale_factors_for_lambda(l, args.alchElecLambdaStart, args.alchVdwLambdaEnd)

    if args.lambda0 is not None and args.lambda1 is not None:
        electro_orig_lambda0, lj_orig_lambda0, \
        electro_orig_lambda1, lj_orig_lambda1 = do_fep_scaling(electro_orig, lj_orig, args.lambda0,
                                                               args.lambda1, args.alchElecLambdaStart,
                                                               args.alchVdwLambdaEnd)
        # TODO: If we have patches, do the FEP scaling with that set of energies as well

    # Generate CSV files and accompanying plots
    for datatype, data in all_data:
        prefix = '%s_%s' % (datatype, args.ligand_resname)
        ligand_atoms = u.select_atoms("resname %s" % args.ligand_resname)
        labels = ['%s %d %.3f' % (a.type, a.index - ligand_atoms[0].index + 1, a.charge) for a in ligand_atoms]

        with open('%s.csv' % prefix, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(labels)
            for row in data:
                writer.writerow(row)
        print >>sys.stderr, 'Wrote energy data to %s.csv' % prefix

        generate_plot(datatype, data, args.dont_plot_atomtypes)
