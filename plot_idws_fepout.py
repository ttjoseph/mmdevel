#!/usr/bin/env python
#
# Plot fepout from NAMD with IDWS on.
# Two plots on a letter-sized paper: raw FEP data, and cumulative sum of ∆E vs lambda
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from plot_many_feps import parse_fepout

def main():
    ap = argparse.ArgumentParser(description='Plot FEP output (.fepout) after being processed with deinterleave_idws.py')
    ap.add_argument('-t', '--title', help='Title to draw above graph')
    ap.add_argument('-o', '--output', help='Output file. Default: <prefix>.pdf')
    ap.add_argument('prefix', help='Prefix to post-process .fepout files. Example: foo -> foo_fwd.fepout and foo_bwd.fepout')
    args = ap.parse_args()

    fwd_fname, bwd_fname = f"{args.prefix}_fwd.fepout", f"{args.prefix}_bwd.fepout"

    fwd, fwd_dg, fwd_lambdas = parse_fepout(fwd_fname)
    bwd, bwd_dg, bwd_lambdas = parse_fepout(bwd_fname)

    fig, ax = plt.subplots(2, 1, figsize=(7.5, 10))
    fig.tight_layout(h_pad=2.0, w_pad=0.2, pad=5)

    ax[0].plot(np.arange(0, len(fwd))/len(fwd), fwd, linewidth=0.2, label='Forward')
    ax[0].plot(np.arange(0, len(bwd))/len(fwd), bwd, linewidth=0.2, label='Backward')
    ax[1].plot(fwd_lambdas, np.cumsum(fwd_dg), marker='.', label='Forward: {:.2f} kcal/mol'.format(np.sum(fwd_dg)))
    ax[1].plot(bwd_lambdas, np.cumsum(bwd_dg), marker='.', label='Backward: {:.2f} kcal/mol'.format(np.sum(bwd_dg)))
    ax[0].set_xlabel(u'Fraction of all timesteps (\u03bb: 0 → 1)')
    ax[0].set_ylabel(u'\u0394\u0394G (kcal/mol)')
    ax[1].set_xlabel(u'\u03bb')
    ax[1].set_ylabel(u'\u0394G (kcal/mol)')

    if args.title is not None:
        ax[0].set_title(args.title, fontsize=8, pad=3)

    out_fname = args.output or f"{args.prefix}.pdf"
    pdf = PdfPages(out_fname)
    pdf.savefig(fig, bbox_extra_artists=[ax[1].legend()])
    pdf.close()
    print(f"Wrote graphs to {out_fname}.")

if __name__ == "__main__":
    main()
