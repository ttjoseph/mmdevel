#!/usr/bin/env python3
import argparse
import sys
import os
from glob import glob
from collections import defaultdict
import numpy as np
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def pretty(s):
    s = s.replace('MDv2', '')
    s = s.replace('MD', '')
    s = s.replace('_try2', '')
    s = s.replace('_run2', '')
    s = s.replace('LigandRMSD/', '')
    s = s.replace('r1', '').replace('r2', '')
    s = s.replace('SKE', 'S-ketamine')
    s = s.replace('SKP', 'S-ketamine (+1)')
    s = s.replace('SNKE', 'S-norketamine')
    s = s.replace('SNKP', 'S-norketamine (+1)')
    s = s.replace('RNKE', 'R-norketamine')
    s = s.replace('RNKP', 'R-norketamine (+1)')
    s = s.replace('TRKE', '(2R,6R)-hydroxynorketamine')
    s = s.replace('TRKP', '(2R,6R)-hydroxynorketamine (+1)')
    s = s.replace('RKE', 'R-ketamine')
    s = s.replace('RKP', 'R-ketamine (+1)')
    s = s.replace('Mu', 'M')
    s = s.replace('Kappa', 'K')
    s = s.replace('.WithDISU', '')
    s = s.replace('HSP297', 'H297+')
    s = s.replace('HSP291', 'H291+')
    s = s.replace('_', ' ')
    # Get rid of extra whitespace
    s = ' '.join(s.split())
    return s



def parse_fepout(fnames, verbose=False):
    fepenergy, deltas, lambdas = list(), list(), list()
    lines = list()

    if isinstance(fnames, str):
        fnames = [fnames,]

    for fname in fnames:
        if verbose:
            print(f"parse_fepout: Processing {fname}...", file=sys.stderr)
        if os.path.exists(fname) is False:
            print(f"parse_fepout: I was asked to parse file {fname} but it doesn't seem to exist",
                file=sys.stderr)
            return None, None, None
        if os.path.getsize(fname) == 0:
            print(f"parse_fepout: File {fname} exists but has zero size", file=sys.stderr)
            return None, None, None        
        with open(fname) as f:
            theselines = f.readlines()
            goodlines = []
            line_counter = 1
            for line in theselines:
                # Check for malformed FepEnergy: lines
                if line.startswith('FepEnergy'):
                    tokens = line.split()
                    if len(tokens) < 10:
                        print(f"parse_fepout: Messed up FepEnergy line in {fname}:{line_counter} with line length {len(line)}", file=sys.stderr)
                        print(f'Here are the first 200 characters in the line ({len(tokens)} tokens):', file=sys.stderr)
                        print(line[:200], file=sys.stderr)
                        continue
                goodlines.append(line)

            lines.extend(goodlines)

    for line in lines:
        tokens = line.split()
        if line.startswith('FepEnergy'):
            # A single sample
            # The energy is always in the tenth column
            if len(tokens) < 10:
                print(f"parse_fepout: Skipping messed up FepEnergy line", file=sys.stderr)
                continue
            energy = float(tokens[9])
            if energy > 20:
                print(f"parse_fepout: Skipping messed up FepEnergy energy {tokens[9]}", file=sys.stderr)
                continue
            fepenergy.append(energy)
        elif line.startswith('#Free'):
            # delta-G at end of window, in the twelfth column
            deltas.append(float(tokens[11]))
        elif line.startswith('#NEW FEP WINDOW'):
            # Lambda value in seventh column
            lambdas.append(float(tokens[6]))

    if len(lambdas) != len(deltas):
        print(f"parse_fepout: While processing {fnames}:")
        print(f"parse_fepout: We have {len(lambdas)} lambdas but {len(deltas)} delta-G values", file=sys.stderr)
        return None, None, None

    return fepenergy, deltas, lambdas


# Return the first matching glob of files
# This will allow us to use (for example) either dis5A000_fwd.fepout, or dis5A???.fepout
def try_globs(*globspecs):
    for globspec in globspecs:
        fnames = glob(globspec)
        if len(fnames) > 0:
            break
    return sorted(fnames)


def find_fepouts(dirnames, fepdirname, prefixes, verbose=False):
    good_dirnames = set()
    fwd_fnames, bwd_fnames = defaultdict(list), defaultdict(list)
    for dirname in dirnames:
        if os.path.exists(dirname) is False or os.path.isdir(dirname) is False:
            print('Could not find directory {}.'.format(dirname), file=sys.stderr)
            continue

        # Try each prefix (dis5A, dis5B, ...)
        for prefix in prefixes.split(','):
            fwd_fname = try_globs(f'{dirname}/{fepdirname}/{prefix}???_fwd.fepout',
                f'{dirname}/{fepdirname}/{prefix}???.fepout')
            bwd_fname = try_globs(f'{dirname}/{fepdirname}/{prefix}???_bwd.fepout')

            # For each directory, there are two dicts: one for forward fepouts, and another for backward
            # fepouts. We require at least a forward fepout. If no backward fepouts, just save an empty list.
            # This ensures that both dicts have the same number of lists of fepouts, so we can iterate
            # through them later properly and not lose which backward fepout sets correspond to which
            # forward fepout sets.
            if len(fwd_fname) == 0:
                # print('Could not find forward fepouts in {}'.format(dirname))
                continue
            else:
                if verbose:
                    print(f"Found forward fepouts with prefix {prefix} in {dirname}", file=sys.stderr)
                fwd_fnames[dirname].append(fwd_fname)
            if len(bwd_fname) == 0:
                if verbose:
                    print('Could not find backward fepouts in {}'.format(dirname), file=sys.stderr)
                bwd_fnames[dirname].append([])
            else:
                bwd_fnames[dirname].append(bwd_fname)

            # Only care about this directory if we found a fepout
            good_dirnames.add(dirname)

    return list(sorted(good_dirnames)), fwd_fnames, bwd_fnames


def main():
    ap = argparse.ArgumentParser(description='Make a big plot of FEP curves')
    ap.add_argument('dirname', nargs='+', help='Directories of simulation systems, with _fwd and _bwd fepouts in namd subdirectories')
    ap.add_argument('--fepdirname', '-f', default='namd', help='Subdirectory under <dirname> containing .fepout files')
    ap.add_argument('--prefix', '-p', default='dis5A,dis5B,dis5C,idws5A', help='fepout file prefixes to try, comma-separated')
    ap.add_argument('--output', '-o', help='Output file')
    args = ap.parse_args()

    good_dirnames, fwd_fnames, bwd_fnames = find_fepouts(args.dirname, args.fepdirname, args.prefix)

    num_feps = len(good_dirnames)
    print('Processing {} FEPs.'.format(num_feps), file=sys.stderr)

    if len(good_dirnames) == 0:
        exit(1)

    rows_per_page = 4
    output_filename = args.output or f"{args.prefix}.fep.pdf"
    print(f"Writing plots to: {output_filename}", file=sys.stderr)
    pdf = PdfPages(output_filename)

    # Ensure there are enough subplots for all the FEPs
    num_cols = 2
    num_rows = num_feps

    num_pages = int(num_rows / rows_per_page)
    if num_rows % rows_per_page > 0:
        num_pages += 1

    fig, ax = None, None
    # mpl.rcParams['xtick.labelsize'] = 6
    plt.rcParams.update({'font.size': 6})
    legends = []

    for i in range(num_feps): # range(6):  # (for debug)
        page = int(i / rows_per_page)
        offset = i % rows_per_page

        if offset == 0:
            print('Starting page {}'.format(page+1))
            fig, ax = plt.subplots(rows_per_page, num_cols, figsize=(10, 7.5))
            fig.tight_layout(h_pad=2.0, w_pad=0.2, pad=5)

        dirname = good_dirnames[i]
        for fepout_i in range(len(fwd_fnames[dirname])):
            fwd, fwd_dg, lambdas = parse_fepout(fwd_fnames[dirname][fepout_i])
            # Some of these files will be corrupt. TODO: Check before this so we don't have empty cells
            if fwd is None:
                continue
            ax[offset, 0].plot(np.arange(0, len(fwd))/len(fwd), fwd, linewidth=0.2, label='Forward')
            #ax[offset][0].set_title(args.prefix, fontsize=8, pad=3)
            # ax[offset][1].set_title('{}/{}'.format(os.getcwd(), dirname), fontsize=8, pad=3)
            ax[offset, 0].set_title(pretty(dirname), fontsize=8, pad=3)
            ax[offset, 1].plot(lambdas, np.cumsum(fwd_dg), marker='.', label='Forward: {:.2f} kcal/mol'.format(np.sum(fwd_dg)))
            # Only have axis labels on the outer edge of the page
            ax[offset, 0].set_xlabel(u'Fraction of all timesteps (\u03bb: 0 → 1)')
            ax[offset, 0].set_ylabel(u'\u0394\u0394G (kcal/mol)')
            ax[offset, 1].set_xlabel(u'\u03bb')
            ax[offset, 1].set_ylabel(u'\u0394G (kcal/mol)')
            ax[offset, 0].label_outer()
            ax[offset, 1].label_outer()
            ax[offset, 1].yaxis.set_label_position("right")
            ax[offset, 1].yaxis.tick_right()

            # dirname should always be in bwd_fnames
            assert(dirname in bwd_fnames)
            if len(bwd_fnames[dirname][fepout_i]) > 0:
                bwd, bwd_dg, lambdas = parse_fepout(bwd_fnames[dirname][fepout_i])
                if bwd is None:
                    continue
                ax[offset][0].plot(np.arange(0, len(bwd))/len(fwd), bwd, linewidth=0.2, label='Backward')
                ax[offset][1].plot(lambdas, np.cumsum(bwd_dg), marker='.', label='Backward: {:.2f} kcal/mol'.format(np.sum(bwd_dg)))

        # Show the legend because the dG sums are in it
        legends.append(ax[offset][1].legend(loc='lower center', bbox_to_anchor=(0.5, 0.95), ncol=2))

        if offset == (rows_per_page - 1) or offset == (num_feps - 1):
            pdf.savefig(fig, bbox_extra_artists=legends)

    pdf.close()

if __name__ == '__main__':
    main()
