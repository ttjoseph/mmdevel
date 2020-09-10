#!/usr/bin/env python
#
# Ensures the list of NAMD fepout files provided as command line arguments span lambda 0 to 1.
import argparse
import sys
from os.path import isfile
from intervaltree import Interval, IntervalTree


def get_fepout_lambda_range(filename):
    # If the file doesn't exist or isn't actually a file, it can't have a lambda
    if isfile(filename) is False:
        return None

    with open(filename) as f:
        lines = f.readlines()
        #Free energy change for lambda window [ 0 0.00137 ] is 0.382175 ; net change until now is 0.382175
        for i in range(len(lines)-1, 0, -1):
            if lines[i].startswith('#Free energy change for lambda window'):
                # print('get_fepout_lambda_range: found it', file=sys.stderr)
                tokens = lines[i].split()
                try:
                    return float(tokens[7]), float(tokens[8])
                except ValueError:
                    print(f"Unable to parse lambda range from the following line in {filename}:")
                    print(lines[i])
                    return None


def main():
    ap = argparse.ArgumentParser(description="Returns a list of fepout filenames covering λ = [0, 1] with sanity checks")
    ap.add_argument('-w', '--num-width', type=int, default=3, help='Width of number field (e.g. dis5A000 would be 3)')
    ap.add_argument('arglist', nargs='+', metavar='arg', help="<prefix> <beg-end> [prefix] [beg-end] ...")
    args = ap.parse_args()

    if len(args.arglist) % 2 != 0:
        print(f'You provided {len(args.arglist)} argument(s) but only an even number of arguments is correct.', file=sys.stderr)
        print('You want to do: <prefix> <beg-end> [prefix] [beg-end] ...')
        exit(1)

    # Parse which fepout files the user wants us to look at
    num_args = len(args.arglist)
    prefixes = args.arglist[0:num_args:2]
    lambda_ranges = args.arglist[1:num_args:2]
    lambda_begin, lambda_end = [], []
    for r in lambda_ranges:
        begin, end = (int(x) for x in r.split('-'))
        lambda_begin.append(begin)
        lambda_end.append(end)

    # Find all the fepout files to examine
    # Extract the lambda interval from each fepout file and add it to the interval tree
    # Along the way, ensure that that each fepout file:
    #   - exists to begin with
    #   - is complete, having a correct last line
    #   - TODO: doesn't have some crazy high-magnitude ∆A at the end
    files = dict()
    lambda_coverage = IntervalTree()
    last_idws_fname, last_idws_window = None, None
    for i in range(len(prefixes)):
        prefix = prefixes[i]
        for num in range(lambda_begin[i], lambda_end[i]+1):
            num_str = f"{num}".zfill(args.num_width)
            fname = f'{prefix}{num_str}.fepout'

            if fname in files:
                print(f"Looks like you specified {fname} twice. I'll allow it.")
            if isfile(fname) is False:
                print(f'{fname} is not a file that is here', file=sys.stderr)
                continue
            lambdas = get_fepout_lambda_range(fname)
            if lambdas is None:
                print(f'{fname} is not a well formed fepout file', file=sys.stderr)
                continue
            files[fname] = lambdas
            if lambdas[0] == 1.0 and lambdas[1] < lambdas[0]:
                last_idws_fname = fname
                last_idws_window = lambdas
            else:
                # print(f'{fname}: {lambdas}', file=sys.stderr)
                if lambda_coverage.overlaps(lambdas[0], lambdas[1]):
                    print(f"Lambdas {lambdas} already seen. Perhaps you are sloppily trying to merge two runs.")
                    return None
                lambda_coverage.addi(lambdas[0], lambdas[1], fname)
    
    # Ensure complete coverage from lambda = 0 to 1
    missing_range = IntervalTree()
    missing_range.addi(0.0, 1.0)
    for interval in lambda_coverage:
        missing_range.chop(interval.begin, interval.end)

    if len(missing_range) > 0:
        print('Error: You are missing some lambdas:', file=sys.stderr)
        print('\n'.join([f"    {x.begin} to {x.end}" for x in sorted(missing_range)]), file=sys.stderr)
        exit(1)


    files_to_print = [x.data for x in sorted(lambda_coverage)]
    if last_idws_window is not None and lambda_coverage.overlap(last_idws_window[1], last_idws_window[0]):
        files_to_print.append(last_idws_fname)
    else:
        print(f"Last (backwards) IDWS window {last_idws_fname} spanning {last_idws_window} is messed up in some way",
            file=sys.stderr)
        exit(1)

    print(' '.join(files_to_print))

if __name__ == "__main__":
    main()