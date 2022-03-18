#!/usr/bin/env python
#
# Ensures the list of NAMD fepout files provided as command line arguments span lambda 0 to 1.
import argparse
import sys
from os.path import isfile
from intervaltree import Interval, IntervalTree
from tqdm import tqdm

def get_fepout_lambda_range(filename):
    # If the file doesn't exist or isn't actually a file, it can't have a lambda
    if isfile(filename) is False:
        return []

    windows = []

    with open(filename) as f:
        lines = f.readlines()
        #Free energy change for lambda window [ 0 0.00137 ] is 0.382175 ; net change until now is 0.382175
        for i in range(len(lines)-1, 0, -1):
            if lines[i].startswith('#Free energy change for lambda window'):
                # print('get_fepout_lambda_range: found it', file=sys.stderr)
                tokens = lines[i].split()
                try:
                    windows.append((float(tokens[7]), float(tokens[8])))
                except ValueError:
                    print(f"Unable to parse lambda range from the following line in {filename}:")
                    print(lines[i])
                    # return None

    return windows


def validate_fepouts(fileslist):
    """Validate a set of fepout files, returning the list back if it's good, and None if not.

    Extract the lambda interval from each fepout file and add it to the interval tree
    Along the way, ensure that that each fepout file:
      - exists to begin with
      - is complete, having a correct last line
      - TODO: doesn't have some crazy high-magnitude ∆A at the end"""
    
    lambda_coverage = IntervalTree()
    last_idws_fname, last_idws_window = None, None
    for fname in tqdm(sorted(list(set(fileslist))), unit=' files'):
        if isfile(fname) is False:
            print(f'{fname} is not a file that is here', file=sys.stderr)
            continue
        windows = get_fepout_lambda_range(fname)
        if windows == []:
            # print(f'{fname} is not a well formed fepout file', file=sys.stderr)
            continue
        for lambdas in windows:
            if lambdas[0] == 1.0 and lambdas[1] < lambdas[0]:
                last_idws_fname = fname
                last_idws_window = lambdas
            else:
                # print(f'{fname}: {lambdas}', file=sys.stderr)
                if lambda_coverage.overlaps(lambdas[0], lambdas[1]):
                    print(f"validate_fepouts: Info: {fname} ({lambdas[0]} to {lambdas[1]}) overlaps with:")
                    for iv in sorted(lambda_coverage.overlap(lambdas[0], lambdas[1]), key=lambda x: x.data):
                        print(f"    {iv.data}: {iv.begin} to {iv.end}")
                    # print(f"Perhaps you are sloppily trying to merge two runs.")
                    # return None
                lambda_coverage.addi(lambdas[0], lambdas[1], fname)

    # Ensure complete coverage from lambda = 0 to 1
    missing_range = IntervalTree()
    missing_range.addi(0.0, 1.0)
    for interval in lambda_coverage:
        missing_range.chop(interval.begin, interval.end)

    if len(missing_range) > 0:
        print('validate_fepouts: Error: You are missing some lambdas:', file=sys.stderr)
        print('\n'.join([f"    {x.begin} to {x.end}" for x in sorted(missing_range)]), file=sys.stderr)
        missing_total = sum([x.end - x.begin for x in missing_range])
        print(f'The calculation appears to be roughly {100*(1 - missing_total):.0f}% done.', file=sys.stderr)
        return None

    files_to_print = [x.data for x in sorted(lambda_coverage)]
    if last_idws_window is not None and lambda_coverage.overlap(last_idws_window[1], last_idws_window[0]):
        files_to_print.append(last_idws_fname)
    else:
        print(f"validate_fepouts: Last (backwards) IDWS window {last_idws_fname} spanning {last_idws_window} is messed up in some way",
            file=sys.stderr)
        return None

    print('validate_fepouts: Complete set of IDWS fepout files is present. Hooray!', file=sys.stderr)
    return(list(sorted(list(set(files_to_print)))))


def main():
    ap = argparse.ArgumentParser(description="Returns a list of IDWS fepout filenames covering λ = [0, 1] with sanity checks")
    ap.add_argument('fileslist', nargs='+')
    args = ap.parse_args()
    files_to_print = validate_fepouts(args.fileslist)
    if files_to_print is not None:
        print(' '.join(files_to_print))


if __name__ == "__main__":
    main()
