#!/usr/bin/env python2.7

# Project "Unravel Motifs" - Align sequences
# DD2404 - NADA, 2015
# Authors: Andreas Jangmo, Johan Lord

from Bio import SeqIO
import argparse, sys, time

# Usage: ./cutseq [-n nucleotides -p #position] [-s #start -e #end] [-l #length] [-o outfile] [fastafile(s)]
# Cutseq takes nucleotide or amino acid sequences in a FASTA formatted
# file and cuts the sequnces between #start and #end positions OR a sequence
# of #length starting from the first character. Option -l can be used
# in conjuction with #start, but #length is ignoredif both #start and #end.
# are given Option -n replaces nucleotides XXX at #position with "-" if present;
# eg if -n is "ATG" it becomes "---" but "ACG" at the same positions remains unaltered.

# If fasta file is omitted, input is read from stdin. Output defaults to
# stdout.

# Example
# You have a FASTA file with a number of sequences aligned so that the 
# start codon "ATG" is at position 16 and you wish to prepare output for creating
# a sequence logo using 15 nucleotides before the start codon and 15 after.
# The file tests/homo_sapiens_ch1_genes_upstream15.fa contains 42 genes from
# homo sapiens chromosome 1.
# Executing...
# ./cutseq.py -n ATG -p 16 -l 33 -r ../tests/homo_sapiens_ch1_genes_upstream15.fa
# OR
# ./cutseq.py -n ATG -p 16 -s 1 -e 33 -r ../tests/homo_sapiens_ch1_genes_upstream15.fa
# gives you a file like this:
# CGAGCGGCCGCCAAC---CTCTTTGAGGGCTTG
# TTTCCGACGGAGTGA---GCGGCGGCGGCTGGG
# ...

class CutSeq:
	def __init__(self, fileResources, format='fasta'):
		"""
		Constructor for the CutSeq class where a file resource or a list
		of file resources containing the nucleotide sequences in specified format.
		
		Arguments:
			fileResources	- a file resource or a list of several file resources
			format			- the format of the nucleotide sequences in the
							file resources
		"""
		if isinstance(fileResources, list) == False: self.fileResources = [fileResources]
		else: self.fileResources = fileResources
		self.length = False
		self.skipSeq = False
		self.start = 0
		self.end = False
		self.n = 0
	
	# Option -l
	def setLength(self, length):
		"""
		Specify length of nucleoutide sequence to extract.
		"""
		self.length = self._toInteger(length, 'length must be integer greater than 0')
		return self
	
	# Option -n and -p
	def setSkip(self, seq, position):
		"""
		Specify if a nucleotide sequence at a certain position is to be omitted, such
		as the start codon "ATG".
		
		Arguments:
			seq		 - nucleotide sequence
			position - position of seq, first character is 1
		"""
		self.skipSeq = seq
		self.skipLen = len(seq)
		self.position = self._toInteger(position, 'position must be integer greater than 0') - 1
		return self
	
	# Option -s and -e
	def setStartEnd(self, start, end=False):
		"""
		Specify starting and ending position of extraction, specifying both overrides
		length if specified. First character is 1.
		
		Arguments:
			start	- position of first character to extract
			end		- position of last character to extract
		"""
		self.start = self._toInteger(start, 'start must be integer greater than 0') - 1
		if end != False:
			self.end = self._toInteger(end, 'end must be integer greater than 0') - 1
			self.length = self.end - self.start + 1
		return self
	
	# Option -d: Is sequence headers to be included?
	def setHeaders(self, state):
		self.headers = state
		return self
	
	# Option -r: Should a report be generated?
	def setReport(self, state):
		self.report = state
		self.startTime = time.time()
		return self
	
	# Option -f: What format is contained in file resources?
	def setFormat(self, format):
		self.format = format
		return self
	
	# Convert 'integer' from string to integer and raise 'error' if conversion
	# not possible or integer is negative
	def _toInteger(self, integer, error):
		try:
			i = int(integer)
			if i <= 0: raise ValueError(error)
			return i
		except Exception:
			raise ValueError(error)
	
	# Display progress, counting number of sequences included
	def _progress(self):
		if self.report:
			self.n += 1
			n = len(str(self.n))
			if self.n > 1: print '\r' * n,
			print self.n,
	
	# This is the workhorse
	def _cut(self, seq):
		seqLen = len(seq)
		if self.length != False: end = self.start + self.length
		else: end = seqLen
		# If the sequence length is shorter that the desired, the sequence
		# is to be skipped.
		if end > seqLen:
			return False
			exit
		# If a short sequence of nucletides is to be skipped, replace these
		# with '-'
		if self.skipSeq != False:
			a = seq[self.start:self.position]
			if seq[self.position:self.position+self.skipLen] == self.skipSeq:
				a += '-' * self.skipLen + seq[self.position+self.skipLen:end]
			else: a += seq[self.position:end]
			return a
		else:
			return seq[self.start:end]
	
	# When object is treated as a string, parse all files, cut all
	# sequences and return the result as a string
	def __str__(self):
		output = ''
		for file in self.fileResources:
			for acc in SeqIO.parse(file, self.format):
				self._progress()
				piece = self._cut(str(acc.seq))
				if self.headers: output += '>' + acc.id + '\n'
				if piece != False: output += piece + '\n'
		if self.report:
			print 'sequences parsed in ' + str(round(time.time()-self.startTime,2)) + ' seconds.'
		return output
		
def main():
	# Set program description
	parser = argparse.ArgumentParser(description="CutSeq cuts out portions from nucleotide or amino acid sequences from FASTA formatted files. 2015 (c) Lord Jangmo Software Curp.")
	# Set options
	parser.add_argument('-r', help='Display process report.', action='store_const', const=True, default=False)
	parser.add_argument('-l', help='Set length of sequence to extract')
	parser.add_argument('-f', help='Set format of input data, can be any accepted by BioPython, defaults to FASTA', default='fasta')
	parser.add_argument('-n', help='Set letters to skip')
	parser.add_argument('-p', help='Start position of letters to skip')
	parser.add_argument('-s', help='Set starting position of extraction')
	parser.add_argument('-e', help='Set ending position of extraction', default=False)
	parser.add_argument('-d', help='Include description headers of sequences', action='store_const', const=True, default=False)
	parser.add_argument('-o', type=argparse.FileType('w'), default=False, help='Direct output to file')
	# Set optional argument infile from which input is read
	parser.add_argument('infile', type=argparse.FileType('r'), nargs='*', default=sys.stdin,
		help='File(s) containing FASTA formatted nucleotide sequences. If omitted input is read from stdin')
	# Parse command line arguments
	args = parser.parse_args()
	try:
		# Create a CutSeq object and set some options
		c = CutSeq(args.infile).setHeaders(args.d).setReport(args.r).setFormat(args.f)
		if args.l != None:
			c.setLength(args.l)
		if args.n != None and args.p != None:
			c.setSkip(args.n, args.p)
		elif (args.n != None and args.p == None) or (args.n == None and args.p != None):
			raise Exception('Both -s and -p need to be specified')
		if args.s != None:
			c.setStartEnd(args.s, args.e)
		# Direct output to file if -o is specified, otherwise stdout
		if args.o != False: args.o.write(str(c))
		else: print c,
	except KeyboardInterrupt:
		print("\nInterrupted by user.\nGoodbye!")
	except Exception as e:
		if str(e).find("urlopen error") != -1: print("Possibly a connection related error. Stupid donkey.")
		else: print(e)

if __name__ == '__main__':
	main()
#print 'a\r',
#print 'b'
#align(file, 'ATG', 15)
