#!/usr/bin/env python

# unravel_motifs, Intron filter
# DD2404 - NADA, 2015
# Authors: Andreas Jangmo, Johan Lord

import re, sys
from Bio import SeqIO

# Usage: ./intronfilter.py genes.fa exons.fa
# where genes.fa is a fasta file of genes with only geneid as header
# and exons.fa is a fasta file of exons with only geneid as header.

# As of 151223/16:55 this code is under construction. What we are
# trying to achieve is to input a fasta file of genes (or transcripts?)
# and get starting and end positions of all exons in each gene so that
# we can easily filter out the area we are interested in (i.e 15 bp upstream
# and 15 bp downstream the beginning and end of the first intron).

# Now the code just prints out the start and end. The problem right now is
# that not all exons are found within the gene.

def intronfilter(gene, exons):
	"""
	As of now, this function prints start and end of each exon
	found in the gene. If the exon is not found it prints 'Exon
	not found'. Preferably all exons should be found so a 'Exon
	not found' signals that something is wrong.

	Arguments:
	  - gene    - a sequence of the entire gene as a string.
	  - exons   - a list of all exons for that gene as strings.
	return: None
	"""
	for exon in exons:
			ptrn = re.compile(exon)
			m = ptrn.match(gene)
			if m:
				print(m.start(), m.end())
			else:
				print('Exon not found!')

def main():
	genes = list(SeqIO.parse(sys.argv[1], "fasta"))
	exons = list(SeqIO.parse(sys.argv[2], "fasta"))
	for gene in genes:
		gene_introns = []
		gene_exons = []
		for exon in exons:
			if gene.id == exon.id:
				gene_exons.append(str(exon.seq))
		print("----GENE: "+gene.id+"----")
		gene_introns = intronfilter(str(gene.seq),gene_exons)

if __name__ == "__main__":
	main()