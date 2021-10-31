import argparse
from collections import defaultdict
from itertools import product
from os.path import realpath, basename
import sys

# Atom types which we shouldn't consider equivalent, even if they have the same L-J parameters
ATOMTYPES_NOT_EQUIVALENT = [
    ('NH3', 'NH2'),
    ('NH3', 'NC2'),
    ('NH3', 'NG2D1'),
    ('NH3', 'NG2O1'),
    ('NH3', 'NG2P1'),
    ('NG3P0', 'NG2D1'),
    ('NG3P0', 'NG2O1'),
    ('NG2P1', 'NG3P0'),
    ('NG2P3', 'NG3P0'),
    ('NG2P3', 'NG3P1'),
    ('NC2', 'NG3P0'),
    ('NG2R43', 'NG3P0'),
    ('NG2R43', 'NH3'),
    ('NG2R50', 'NG3P0'),
    ('NG2R50', 'NH3'),
    ('NG2R51', 'NG3P0'),
    ('NG2R51', 'NG2P1'),
    ('NG2R61', 'NG3P0'),
    ]

ATOMTYPES_EQUIVALENT = [
    ('CTL1', 'CG331'),
    ('CTL2', 'CG321'),
    ('CTL3', 'CG331'),
    ('CTL5', 'CG334'),
    ('CA', 'CG2R61'),
    ('CPT', 'CG2RC0'),
    ('NC2', 'NG2P1'),
    ('NTL', 'NG3P0'),
    ('NH3L', 'NH3'),
    ('NH3L', 'NG3P3'),
    ('NH3', 'NG3P3'),
    ('PL', 'PG1'),
]


def keep_only_equivalent(atomtype, equiv_atomtypes):
    """Returns equiv_atomtypes but only those atom types which are equivalent to
    atomtype, as determined by ATOMTYPES_EQUIVALENT.

    Perhaps there is a more efficient way to accomplish the overall purpose.
    """
    global ATOMTYPES_EQUIVALENT
    out_equiv_atomtypes = list()

    for a, b in ATOMTYPES_EQUIVALENT:
        if atomtype == a and b in equiv_atomtypes:
            print(f"Kept ({a}, {b})", file=sys.stderr)
            out_equiv_atomtypes.append(b)
        elif atomtype == b and a in equiv_atomtypes:
            print(f"Kept ({b}, {a})", file=sys.stderr)
            out_equiv_atomtypes.append(a)
    # An atomtype is equivalent to itself
    out_equiv_atomtypes.append(atomtype)
    return out_equiv_atomtypes


def cull_not_equivalent(atomtype, equiv_atomtypes):
    """Returns equiv_atomtypes modulo those atom types which are not equivalent to atomtype,
    as determined by ATOMTYPES_NOT_EQUIVALENT.
    """
    global ATOMTYPES_NOT_EQUIVALENT

    out_equiv_atomtypes = list(set(equiv_atomtypes))

    for a, b in ATOMTYPES_NOT_EQUIVALENT:
        if atomtype == a and b in out_equiv_atomtypes:
            print(f"Culled {a} due to {b}", file=sys.stderr)
            out_equiv_atomtypes.remove(b)
        elif atomtype == b and a in out_equiv_atomtypes:
            print(f"Culled {b} due to {a}", file=sys.stderr)
            out_equiv_atomtypes.remove(a)
    
    return out_equiv_atomtypes


def get_lj_from_prm(fname, ljparams=defaultdict(list), atomtypes={}, from_fname={}):
    """Extracts atom type names and parameter values from CHARMM parameter/stream files."""

    section_types = ['ATOM', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'NONB', 'NBFI']
    current_section = None
    with open(fname) as f:
        for line in f.readlines():
            line = line.strip()
            # Strip comments and tokenize
            if '!' in line:
                line = line[:line.index('!')].strip()
            l = line.split()
            # Ignore blank lines
            if len(l) == 0:
                continue
            if l[0].lower() == 'cutnb':
                continue
            # Decide what section we're in (ATOM, BOND, NONB, ...)
            is_new_section = False
            for section_type in section_types:
                if l[0].startswith(section_type):
                    current_section = section_type
                    is_new_section = True
            # Don't treat the section header as a line with data in it
            if is_new_section:
                continue

            if current_section == 'NONB':
                # Remember atom type and its L-J parameters
                atomtype, epsilon, rmin2, epsilon14, rmin214 = None, None, None, None, None
                if len(l) >= 4:
                    atomtype, _, epsilon, rmin2 = l[:4]
                if len(l) >= 7:
                    _, _, _, _, _, epsilon14, rmin214 = l[:7]
                # We parse into floats because '1.85' != '1.8500'
                try:
                    epsilon = float(epsilon) if epsilon is not None else ''
                    rmin2 = float(rmin2) if rmin2 is not None else ''
                except ValueError:
                    # If we can't even parse the first two parameters into floats,
                    # give up on this line
                    print('Could not make sense of this (which might be OK):', line, file=sys.stderr)
                    continue
                epsilon14 = float(epsilon14) if epsilon14 is not None else ''
                rmin214 = float(rmin214) if rmin214 is not None else ''

                p = f"{epsilon},{rmin2},{epsilon14},{rmin214}"
                atomtypes[atomtype] = p
                ljparams[p].append(atomtype)
            elif current_section == 'NBFI' and len(l) >= 4:
                # Handle NBFIX entries
                # We remember these so our new NBFIX entries don't step on them
                atomtype1, atomtype2, epsilon, rmin2 = l[:4]
                p = f"{epsilon},{rmin2}"
                atomtype = f"{atomtype1},{atomtype2}"
                atomtypes[atomtype] = p
                ljparams[p].append(atomtype)
                from_fname[atomtype] = fname

    return ljparams, atomtypes, from_fname


def main():
    ap = argparse.ArgumentParser(description='Generate new NBFIXes with equivalent atom types')
    ap.add_argument('target')
    ap.add_argument('template', nargs='+')
    args = ap.parse_args()

    # Extract L-J atom parameters from all these files
    ljparams, atomtypes, from_fname = defaultdict(list), {}, {}
    for fname in args.template:
        if realpath(fname) == realpath(args.target):
            print('Target filename is included in template filename list. You don\'t want that.',
                file=sys.stderr)
            exit(1)
        ljparams, atomtypes, from_fname = get_lj_from_prm(fname, ljparams, atomtypes, from_fname)


    # NBFIX lines have the following form:
    # CG2R51 CG334    -0.04205    4.3150
    # We want to duplicate the same NBFIX for all L-J-equivalent pairs

    # Go through all NBFIX atom types and find all equivalent atom types from template .prm files
    ljparams_nbfix, atomtypes_nbfix, _ = get_lj_from_prm(args.target)
    # print(atomtypes.keys())

    template_fnames_str = '\n'.join([f"*   {basename(fname)}" for fname in args.template])
    not_equiv_str = '\n'.join([f"*   {a} != {b}" for a, b in ATOMTYPES_NOT_EQUIVALENT])
    is_equiv_str = '\n'.join([f"*   {a} == {b}" for a, b in ATOMTYPES_EQUIVALENT])
    
    print(f"""* NBFIX entries generated from: {basename(args.target)}
* Ground rules:
{is_equiv_str}
* Force field files used:
{template_fnames_str}

NBFIX
! atomtype1 atomtype2 epsmin rmin""")
    lines = []
    for k, v in atomtypes_nbfix.items():
        if k in atomtypes: # also check for a2,a1
            print(f"Atom types {k} already present in some NBFIX in the template parameter files.",
                file=sys.stderr)
            exit(1)
        a1, a2 = k.split(',')
        epsilon, rmin = v.split(',')
        missing_atomtype = None
        if a1 not in atomtypes:
            missing_atomtype = a1
        elif a2 not in atomtypes:
            missing_atomtype = a2
        if missing_atomtype is not None:
            print(f"Error: Atom type {missing_atomtype} in {args.target} not present in template parameter files.",
                file=sys.stderr)
            print(f"       You probably haven't specified all the relevant files for the force field of interest.",
                file= sys.stderr)
            exit(1)
        equiv_a1 = ljparams[atomtypes[a1]]
        equiv_a2 = ljparams[atomtypes[a2]]

        equiv_a1 = keep_only_equivalent(a1, equiv_a1)
        equiv_a2 = keep_only_equivalent(a2, equiv_a2)
        # equiv_a1 = cull_not_equivalent(a1, equiv_a1)
        # equiv_a2 = cull_not_equivalent(a2, equiv_a2)

        print(f'{a1}: {equiv_a1}', file=sys.stderr)
        print(f'{a2}: {equiv_a2}', file=sys.stderr)

        # We now have a set of parameters that apply to a bunch of atom type pairs,
        # which can be calculated by cartesian product of lists, so now we repeat those
        # parameters for them all.
        for this_a1, this_a2 in product(equiv_a1, equiv_a2):
            if f'{this_a1},{this_a2}' in atomtypes or f'{this_a2},{this_a1}' in atomtypes:
                if f'{this_a1},{this_a2}' in from_fname:
                    fname = from_fname[f'{this_a1},{this_a2}']
                else:
                    fname = from_fname[f'{this_a2},{this_a1}']
                print(f'Warning: Skipped {this_a1},{this_a2} as there already is an NBFIX for that in {basename(fname)}', file=sys.stderr)
            else:
                idx = f"{a1},{a2}"
                lines.append(f"{this_a1:8s} {this_a2:8s} {epsilon:8s} {rmin:8s} ! {a1} {a2} in {basename(args.target)}")

    # Report multiple entries for the same pair.
    # This can happen when atom types have the same L-J parameters but occur in different NBFIX entries.
    # In this case the user needs to decide which atom types are really equivalent.
    output_lines = []
    prev_atomtype_pair = None
    num_duplicates = 0
    for line in sorted(set(lines)):
        l = line.split()
        if l[:2] == prev_atomtype_pair:
            num_duplicates += 1
        elif num_duplicates > 0:
            output_lines.append(f'! ^^^ {num_duplicates+1} entries for same atom type pair above ^^^')
            num_duplicates = 0
        output_lines.append(line)
        prev_atomtype_pair = l[:2]

    if num_duplicates > 0:
        output_lines.append(f'! ^^^ {num_duplicates+1} entries for same atom type pair above ^^^')
    print('\n'.join(output_lines))
    

if __name__ == '__main__':
    main()
