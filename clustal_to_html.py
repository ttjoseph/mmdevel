#!/usr/bin/env python3
#
# Converts CLUSTAL multiple sequence alignment output to a friendlier HTML output
# because it's a real pain to figure out what specific resid maps to what
import argparse
from Bio import AlignIO


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('alignment', help="CLUSTAL alignment file")
    args = ap.parse_args()

    with open(args.alignment) as f:
        aln = AlignIO.read(f, 'clustal')

    print(f"""<!doctype html>
<html>
<head><title>{args.alignment}</title>
<style>
tr:nth-child(even) {{
  background-color: #f2f2f2;
}}
thead {{
  position: sticky;
  top: 0.5em;
  background: white;
}}
th.rotate {{
  /* Something you can count on */
  height: 140px;
  white-space: nowrap;
}}

th.rotate > div {{
  transform: 
    /* Magic Numbers */
    translate(5px, 51px)
    /* 45 is really 360 - 45 */
    rotate(315deg);
  width: 30px;
}}
th.rotate > div > span {{
  border-bottom: 1px solid #ccc;
  padding: 5px 10px;
}}
table {{ border-spacing: 0px }}
</style>
</head>
<body>
<table>
""")
    print("<thead><tr>")
    for seqid in range(len(aln)):
        print(f'<th scope="col" class="rotate"><div></span>{aln[seqid].id}</span></div></th>')
    print("</tr></thead>")

    idx_count = [0] * len(aln)

    for i in range(len(aln[0])):
        # print(i, aln[0][i], aln[1][i])
        print("<tr>")
        for seqid in range(len(aln)):
            if aln[seqid][i] != '-':
                idx_count[seqid] += 1
                print(f"<td>{aln[seqid][i]}{idx_count[seqid]}</td>")
            else:
                print("<td></td>")
        print("</tr>")


    print(f"""</table></body></html>""")

if __name__ == '__main__':
    main()

