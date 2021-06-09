#!/usr/bin/env python
#
# Plot fepout from NAMD with IDWS on.
# Two plots on a letter-sized paper: raw FEP data, and cumulative sum of ∆E vs lambda
import sys
import argparse
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from libttj import parse_fepout

def main():
    ap = argparse.ArgumentParser(description='Plot FEP output (.fepout) after being processed with deinterleave_idws.py')
    ap.add_argument('-t', '--title', help='Title to draw above graph')
    ap.add_argument('-l', '--labels', help='Comma-separated labels for prefixes, respectively')
    ap.add_argument('-o', '--output', help='Output file. Default: <prefix>.pdf')
    # ap.add_argument('-c', '--cumul-only', type=bool, default=False, help='Plot only cumulative ∆G sum')
    ap.add_argument('prefix', nargs='+', help='Prefix to post-process .fepout files. Example: foo -> foo_fwd.fepout and foo_bwd.fepout')
    args = ap.parse_args()

    labels = args.labels.split(',') if args.labels is not None else None
    fig, ax = plt.subplots(2, 1, figsize=(7.5, 10))

    i = 0
    for prefix in args.prefix:
        fwd_fname, bwd_fname = f"{prefix}_fwd.fepout", f"{prefix}_bwd.fepout"

        # fwd and bwd here are the raw ddG data at each snapshot as recorded
        fwd, fwd_dg, fwd_lambdas = parse_fepout(fwd_fname)
        bwd, bwd_dg, bwd_lambdas = parse_fepout(bwd_fname)

        both_dg = np.concatenate((fwd_dg, bwd_dg), axis=None)
        mean_sum_dg = np.mean([np.abs(np.sum(fwd_dg)), np.abs(np.sum(bwd_dg))])
        hysteresis = np.sum(both_dg)

        this_label = f"{labels[i]}: " if labels is not None else ''
        # Real data cumulative
        ax[0].plot(fwd_lambdas, fwd_dg.cumsum(), alpha=0.5, marker='.', label=f'{prefix} forward')
        ax[0].plot(fwd_lambdas, -bwd_dg.cumsum(), alpha=0.5, marker='.', label=f'{prefix} backward')
        ax[1].plot(fwd_lambdas, np.array(fwd_dg) + np.array(bwd_dg), marker='.', label=f'{prefix} forward + backward')
        data = list(zip(fwd_lambdas, np.array(fwd_dg), np.array(bwd_dg),
            np.array(fwd_dg) + np.array(bwd_dg)))
        for foo in data:
            print(foo)
        fig.text(0.13, 1-0.015*i-0.05, (f"{prefix}: {this_label}mean {mean_sum_dg:.1f} kcal/mol, diff {hysteresis:.1f} kcal/mol"))

        i += 1

    if args.title is not None:
        fig.suptitle(args.title)

    for this_ax in ax:
        this_ax.set_xlabel(u'\u03bb')
        this_ax.set_ylabel(u'\u0394G (kcal/mol)')
        this_ax.set_xticks(np.arange(0, 1.0, 0.1))
        this_ax.legend()
        this_ax.grid()

    ax[0].set_title(f'Cumulative ∆G ({len(fwd_lambdas)} windows)')
    ax[1].set_title('Sum of IDWS forward and backward ∆G values')

    fig.tight_layout(h_pad=2.0, w_pad=0.2, pad=5)
    out_fname = args.output or f"{Path(args.prefix[0]).stem}.pdf"
    pdf = PdfPages(out_fname)
    pdf.savefig(fig)
    pdf.close()
    print(f"Wrote graphs to {out_fname}.")

if __name__ == "__main__":
    main()
 
