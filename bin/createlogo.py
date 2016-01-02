#!/usr/bin/env python

# Project "Unravel Motifs" - Create Logo
# DD2404 - NADA, 2015
# Authors: Andreas Jangmo, Johan Lord

# Usage ./createlogo [-y y axis scale] [-t title for logo] [-o outfile] [infile [infile ...]]

# Example: You have one or several alignment file(s) in flat format (i.e the sequences
# (of the same length) are just listed on each line in a plain text file without names).
# From this you want to create ONE sequence logo named outlogo.png. Lets say the infiles
# are named myinfile_1.txt through myinfile_10.txt and you want y axis to go up to 1.5 and
# the logo to be named "My precious logo", then you would run:
# ./createlogo -y 1.5 -t "My precious logo" -o outlogo.png myinfile_*.txt

import sys, argparse
from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def createlogo(infiles, outfile, title, scale):
	"""
	Creates a sequence logo in PNG format from a textfile of 
	aligned sequences in flat format.

	Arguments:
	  - infile    - name of the flat format alignment file.
	  - outfile   - name of the png logo file.
	  - title     - title for the logo.
	  - scale     - yaxis scale, usually 1.0
	return: None
	"""
	# input a text file with all sequences aligned in flat format.
	for fil in infiles:
		instances = []
		for line in fil:
			# for each line in the file create a seq-object from it
			# and append this object to a list.
			instances.append(Seq(line.strip(), alphabet=IUPAC.ambiguous_dna))
	# create a motif-object. The last sequence in instances
	# is empty, so skip it.
	m = motifs.create(instances[:-1], alphabet=IUPAC.ambiguous_dna)
	# create the sequence logo
	m.weblogo(outfile, logo_title=title, yaxis_scale=scale, stack_width="large")#,alphabet='ambiguous_dna_alphabet')

def main():
	parser = argparse.ArgumentParser(description="createlogo uses the Berkeley weblogo service to create a sequence logo from a flat format text file of aligned sequences.")
	parser.add_argument('-y', help='Set y axis scale.', default = 1.0)
	parser.add_argument('-t', help='Set logo title.', default = "")
	# Set optional argument infile from which input is read
	parser.add_argument('-o', default="outlogo.png", help='Name of the PNG format output file')
	parser.add_argument('infile', type=argparse.FileType('r'), nargs='*', default=sys.stdin,
		help='File(s) in flat format with sequences aligned. If omitted input is read from stdin.')
	# Parse command line arguments
	args = parser.parse_args()
	try:
		createlogo(args.infile, args.o, args.t, args.y)
	except IOError as e2:
		print("Unable to open file: " + e2)
	except KeyboardInterrupt:
		print("\nGoodbye!")
	except Exception as e3:
		print(e3)

if __name__ == '__main__':
	main()