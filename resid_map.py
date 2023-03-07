#!/usr/bin/python
# -*- coding: utf8
from io import StringIO
import argparse
import sys
from Bio import AlignIO

def annotate_numbers(seq):
	ids = []
	counter = 0
	for char in seq:
		if char == '-':
			ids.append(None)
		else:
			ids.append(counter)
			counter += 1
	return ids


def translate(seq_have, seq_want, query_id, noisy=True):
	# Annotate each residue with its original resid, starting from 1
	# Yeah, we do this every time. Sue me.
	id_have = annotate_numbers(seq_have)
	id_want = annotate_numbers(seq_want)

	# The user gave us an original resid, so locate the offset in the array
	# which should give us the equivalent residue in both sequences
	array_idx = id_have.index(query_id)
	if array_idx is None:
		print('# Cannot find that resid in sequence you have', file=sys.stderr)
		return None
	if id_want[array_idx] is None:
		print('# Resid {} does not exist in the other sequence'.format(query_id + 1), file=sys.stderr)
		return None

	if noisy:
		print('{}{} → {}{}'.format(seq_have[array_idx], id_have[array_idx] + 1, seq_want[array_idx], id_want[array_idx] + 1))
	return seq_have[array_idx], id_have[array_idx] + 1, seq_want[array_idx], id_want[array_idx] + 1


def main():
	ap = argparse.ArgumentParser(description='Translate residue IDs according to alignment')
	ap.add_argument('i_have', help='Sequence you have, a or b, defined by their order in the alignment')
	ap.add_argument('resid', nargs='+', type=int, help='Residue ID you have')
	ap.add_argument('--tcl', action='store_true', help='Generate Tcl code that you might use in a VMD script')
	ap.add_argument('--alignment', help="CLUSTAL alignment output")
	args = ap.parse_args()
	seq_a, seq_b = '', ''
	global data

	if args.alignment is not None:
		data = open(args.alignment)
	else:
		data = StringIO(data)

	aln = AlignIO.read(data, 'clustal')

	if args.i_have == 'a':
		seq_have, seq_want = aln[0], aln[1]
	elif args.i_have == 'b':
		seq_have, seq_want = aln[1], aln[0]
	else:
		print('Hey! Specify sequence a or b.', file=sys.stderr)
		exit()

	print('# {} → {}:'.format(seq_have.id, seq_want.id))
	if args.tcl:
		print('array set a_to_b {}')
		print('array set b_to_a {}')

	for resid in args.resid:
		result = translate(seq_have, seq_want, resid - 1, noisy=not args.tcl)
		if args.tcl and result:
			resname1, idx1, resname2, idx2 = result
			print(f"""set a_to_b({idx1}) {idx2}""")
			print(f"""set b_to_a({idx2}) {idx1}""")



# Sample data for ADRA2A homology model
data = """CLUSTAL O(1.2.4) multiple sequence alignment


sp|P08913|ADA2A_HUMAN      MGSLQPDAGNASWNGTEAPGGGARATPYSLQVTLTLVCLAGLLMLLTVFGNVLVIIAVFT	60
homology                   -------------------------------VTLTLVCLAGLLMLLTVFGNVLVIIAVFT	29
                                                          *****************************

sp|P08913|ADA2A_HUMAN      SRALKAPQNLFLVSLASADILVATLVIPFSLANEVMGYWYFGKAWCEIYLALDVLFCTSS	120
homology                   SRALKAPQNLFLVSLASADILVATLVIPFSLANEVMGYWYFGKAWCEIYLALDVLFCTSS	89
                           ************************************************************

sp|P08913|ADA2A_HUMAN      IVHLCAISLDRYWSITQAIEYNLKRTPRRIKAIIITVWVISAVISFPPLISIEKKGGGGG	180
homology                   IV-LCAISLDRYWSITQAIEYNLKRTPRRIKAIIITVWVISAVISFPPLISIEKKGGGGG	148
                           ** *********************************************************

sp|P08913|ADA2A_HUMAN      PQPAEPRCEINDQKWYVISSCIGSFFAPCLIMILVYVRIYQIAKRRTRVPPSRRGPDAVA	240
homology                   PQPAEPRCEINDQKWYVISSCIGSFFAPCLIMILVYVRIYQIAKR---------------	193
                           *********************************************               

sp|P08913|ADA2A_HUMAN      APPGGTERRPNGLGPERSAGPGGAEAEPLPTQLNGAPGEPAPAGPRDTDALDLEESSSSD	300
homology                   ------------------------------------------------------------	193
                                                                                       

sp|P08913|ADA2A_HUMAN      HAERPPGPRRPERGPRGKGKARASQVKPGDSLPRRGPGATGIGTPAAGPGEERVGAAKAS	360
homology                   ------------------------------------------------------------	193
                                                                                       

sp|P08913|ADA2A_HUMAN      RWRGRQNREKRFTFVLAVVIGVFVVCWFPFFFTYTLTAVGCSVPRTLFKFFFWFGYCNSS	420
homology                   --RGRQNREKRFTFVLAVVIGVFVVCWFPFFFTYTLTAVGCSVPRTLFKFFFWFGYCNSS	251
                             **********************************************************

sp|P08913|ADA2A_HUMAN      LNPVIYTIFNHDFRRAFKKILCRGDRKRIV	450
homology                   LNPVIYTIFN-DFRRAFKKILCRGDRKRIV	280
                           ********** *******************
"""

if __name__ == '__main__':
    main()
