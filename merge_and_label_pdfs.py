#!/usr/bin/env python
#
# Merges a bunch of PDF files together, stamping a header of the file name and global page number.
# Useful when you've generated a bunch of individual PDFs of figures and need them all in one file.

import io
import sys
import argparse
import PyPDF2 as PDF
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch


def make_pdf_of_text(text, for_page):
	text_buf = make_pdfbytes_of_text(text, for_page.mediaBox)
	return PDF.PdfFileReader(text_buf).getPage(0)


def make_pdfbytes_of_text(text, media_box):
	buf = io.BytesIO()
	c = canvas.Canvas(buf, pagesize=media_box)
	# Put the text at the top center
	c.drawCentredString((float(media_box.lowerRight[0]) - float(media_box.lowerLeft[0]))/2,
		float(media_box.upperLeft[1]) - 12, text)
	c.save()
	buf.seek(0)
	return buf


def main():
	ap = argparse.ArgumentParser(description="Concatenate PDFs, labeling them with prettified versions of their filenames")
	ap.add_argument('--out', '-o', default='out.pdf', help='Output filename')
	ap.add_argument('--page-number-prefix', '-p', default='', help='Prefix for each page number (e.g. "S" for S1, S2, ...)')
	ap.add_argument('pdfs', nargs='+', help='PDF files to concatenate, in order')
	args = ap.parse_args()

	output = PDF.PdfFileWriter()

	out_page_num = 1
	for fname in args.pdfs:
		# We don't use a "with" block here because PyPDF4 needs each file to remain open until done writing output
		f = open(fname, 'rb')
		this_pdf = PDF.PdfFileReader(f)
		this_num_pages = this_pdf.getNumPages()
		print(f"Working on {fname} ({this_num_pages} pages)...", file=sys.stderr)
		for page_i in range(this_num_pages):
			page = this_pdf.getPage(page_i)

			text_string = f"{fname.replace('_', ' ').replace('.pdf', '')} (p. {args.page_number_prefix}{out_page_num})"
			page.mergePage(make_pdf_of_text(text_string, page))
			output.addPage(page)
			out_page_num += 1

	print(f"Writing output to {args.out}.", file=sys.stderr)
	with open(args.out, 'wb') as f:
		output.write(f)


if __name__ == '__main__':
	main()
