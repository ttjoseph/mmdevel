#!/usr/bin/env python
import sys
import argparse
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import cycler
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from ffparam_charges_vs_qm import load_prm, COULOMB


def calc_mm_interaction_energy(u, ligand_spec, prm, solute_spec='protein or nucleic'):
    """
    Given an MDAnalysis universe, calculate the pairwise interaction energies between the specified
    ligand and the solute.
    :param u: MDAnalysis universe
    :param ligand_spec: Specification to find the ligand molecule, such as "resname APM"
    :param prm: Dict containing Lennard-Jones parameters by atom type, as read from ffparam_charges_vs_qm.load_prm
    :param solute_spec: Specification to find the solute, such as "protein or nucleic"
    :return: A tuple containing electrostatic and Lennard-Jones energies, as NxM arrays, with N ligand atoms and M frames
    """
    ligand_atoms = u.select_atoms(ligand_spec)
    solute_atoms = u.select_atoms(solute_spec)
    ligand_rmin2 = np.array([prm[a.type]['rmin2'] for a in ligand_atoms])
    ligand_epsilon = np.array([prm[a.type]['epsilon'] for a in ligand_atoms])
    solute_rmin2 = np.array([prm[a.type]['rmin2'] for a in solute_atoms])
    solute_epsilon = np.array([prm[a.type]['epsilon'] for a in solute_atoms])

    # Precompute the pairwise charge products, which we will later scale by inverse distance,
    # as well as part of the Lennard-Jones calculation.
    # We can do this otherwise horrible O(n^2) calculation because we know we have a relatively
    # small number of ligand atoms, so the size of the resulting matrix is not cataclysmic.
    # There has got to be a Numpy-array-based way to do this, but since we only do it once
    # we can tolerate the slowness.
    charge_prod = np.zeros((len(ligand_atoms), len(solute_atoms)))
    lj6_partial, eps_partial = np.zeros_like(charge_prod), np.zeros_like(charge_prod)

    for ligand_i in range(len(ligand_atoms)):
        ligand_charge = COULOMB * ligand_atoms[ligand_i].charge
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
        lj6 = lj6_partial / (dists**6)
        frame_lj = np.sum(eps_partial * (lj6**2 - 2*lj6), 1)

        electro.append(frame_electro)
        lj.append(frame_lj)

    return electro, lj


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Calculates single point energies with MM, with respect to a ligand.')
    ap.add_argument('ligand_resname', help='Ligand residue name (e.g. SKE, APM)')
    ap.add_argument('namdconf', help='NAMD configuration containing "parameters" and "structure" keywords from which to glean PSF and force field filenames')
    ap.add_argument('dcd', nargs='+', help='DCD trajectories to go with the PSF file specified in namdconf')
    args = ap.parse_args()

    psf_filename, prm_filenames = None, []

    # Glean PSF and force field parameters from NAMD simulation input file
    print >>sys.stderr, 'Reading NAMD simulation config', args.namdconf
    with open(args.namdconf) as f:
        psf_re = re.compile(r'structure\s+([\w/.]+)')
        prm_re = re.compile(r'parameters\s+([\w/.]+)')
        for line in f:
            m = psf_re.match(line.strip())
            if m: psf_filename = m.group(1)
            m = prm_re.match(line.strip())
            if m: prm_filenames.append(m.group(1))

    print >>sys.stderr, 'Guessed PSF filename:', psf_filename
    print >>sys.stderr, 'Guessed force field parameter files:', ', '.join(prm_filenames)

    u = mda.Universe(psf_filename, args.dcd)

    prm = {}
    for prm_fname in prm_filenames:
        prm.update(load_prm(prm_fname))

    print >>sys.stderr, 'Loaded %d MM parameters' % len(prm)

    ligand_spec = 'resname %s' % args.ligand_resname
    ligand_atoms = u.select_atoms(ligand_spec)
    labels = ['%d %s %.3f' % (a.index - ligand_atoms[0].index + 1, a.type, a.charge) for a in ligand_atoms]

    electro, lj = calc_mm_interaction_energy(u, ligand_spec, prm)

    # Try not to repeat line styles
    pretty_cycler = cycler.cycler(lw=[0.3, 0.8, 1.4]) * cycler.cycler('ls', ['-', ':']) * plt.rcParams['axes.prop_cycle']
    plt.rc('axes', prop_cycle=pretty_cycler)

    for datatype, data in (('electro', electro), ('lj', lj)):
        prefix = '%s_%s' % (datatype, args.ligand_resname)

        with open('%s.csv' % prefix, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(labels)
            for row in data:
                writer.writerow(row)
        print >>sys.stderr, 'Wrote energy data to %s.csv' % prefix

        d = pd.read_csv('%s.csv' % prefix)
        plt.figure(figsize=(8, 4.5)) # 16:9 aspect ratio
        plt.autoscale(tight=True)
        plt.plot(d)
        plt.title('%s energy for %s (%s)' % (datatype.capitalize(), args.ligand_resname, args.namdconf))
        plt.xlabel('Trajectory frame')
        plt.ylabel('Energy (kcal/mol)')
        plt.legend(labels, fontsize=4, loc='upper left', bbox_to_anchor=(1.0, 1.15))
        plt.savefig('%s.pdf' % prefix)
        print >>sys.stderr, 'Wrote a pretty picture to %s.pdf' % prefix
