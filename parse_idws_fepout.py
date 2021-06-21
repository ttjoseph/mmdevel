import sys
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from alchemlyb.parsing.namd import extract_u_nk
from alchemlyb.estimators import BAR
from validate_fepouts import validate_fepouts

GAS_CONSTANT = 1.98720425864083e-3 # kcal⋅K−1⋅mol−1

def main():
    ap = argparse.ArgumentParser(description='Parse fepouts using alchemlyb')
    ap.add_argument('fepout', nargs='+', help='.fepout files, possibly with IDWS energies')
    ap.add_argument('-t', '--temperature', type=float, default=303.15, help='Temperature in K at which simulation was run')
    args = ap.parse_args()

    # TODO: Use validate_fepouts.py algorithm to ensure we have complete set of fepouts
    valid_fepouts = validate_fepouts(args.fepout)
    if valid_fepouts is None:
        exit(1)

    print(f'Loading {len(valid_fepouts)} fepout files...', file=sys.stderr)
    all_data = extract_u_nk(valid_fepouts, args.temperature)
    print(all_data.columns)

    bar_fit = BAR().fit(all_data)
    print('Energies (kcal/mol):')
    print(np.array(bar_fit.delta_f_.iloc[0])*GAS_CONSTANT*args.temperature)


if __name__ == '__main__':
    main()
