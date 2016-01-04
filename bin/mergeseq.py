#!/usr/bin/env python2.7
from Bio import SeqIO
import math,argparse

class Sequences:
	def __init__(self, fileResources, format='fasta'):
		if isinstance(fileResources, list) == False: self.fileResources = [fileResources]
		else: self.fileResources = fileResources
		self.format = format
		self.meta = {}
		self.positions = {}
	
	# Set position in headers
	def setPos(self, name, position):
		self.__dict__[name] = position - 1
		return self
	
	# Iterate through sequences and return a tuple with filename and SeqBio object
	def _parse(self):
		for fr in self.fileResources:
			for seq in SeqIO.parse(fr, self.format):
				yield fr.name,seq
	
	# Only parse headers in the provided sequences
	def parseMeta(self):
		for name,seq in self._parse():
			parts = seq.id.split('|')
			try:
				coord = [(int(parts[self.start]),int(parts[self.end]))]
				id = parts[self.id]
				# Store strand of of sequence
				self.positions[id] = int(parts[self.strand])
				# Add coordinate for current sequence in list
				if id in self.meta: self.meta[id] += coord
				else: self.meta[id] = coord
			except:
				raise Exception('Error parsing header '+seq.id+' in file '+name+'\nDid you specify the header structure properly?')
		if len(self.meta) == 0: raise Exception('No sequences where found in subsequence file!')
		return self
	
	# Create unions and complements according to coordinates stored in self.meta
	def regions(self, gene):
		coordinates = self.meta[gene]
		coordinates.sort()
		unionCoord = []
		complCoord = []
		coordinates.sort()
		start, end = coordinates[0]
		for coord in coordinates[1:]:
			if coord[0] > end:
				unionCoord += [(start,end)]
				complCoord += [(end+1,coord[0])]
				start, end = coord
			if coord[0] < start: start = coord[0]
			if coord[1] > end: end = coord[1]
		unionCoord += [(start,end)]
		return unionCoord,complCoord

class MergeSeq:
	def __init__(self, geneResources, exonResources, format='fasta'):
		self.seqs = Sequences(geneResources, format)
		self.subseqs = Sequences(exonResources, format)
		self.format = format
		self.errors = []
	
	# Pad string with a character to desired width
	def _pad(self, seq, start, geneLen, pad='-'):
		endAdd = geneLen - start - len(seq)
		return pad * start + seq + pad * endAdd + '\n'

	# Convert string to FASTA format
	def _toFASTA(self, seq, width=60):
		o = ''
		end = '\n'
		a = 0
		b = 0
		seqLen = len(seq)
		# r equals the number of rows with 'width' that are to be created
		r = int(math.floor(seqLen/60)) + 1
		# If the sequence fills 'r' rows entirely, no line ending needs to be appended at the return statement
		if (r-1)*60 == seqLen: end = ''
		for x in range(1, r):
			b = 60 * x
			o += seq[a:b] + '\n'
			a = b
		return o + seq[a:] + end
	
	# Option -hs and -hss
	def setStructure(self, seq, subseq):
		"""
		Set the structure of the headers in FASTA file(s). seq(uences) and subseq(ences) needs to have
		position of ID, start and end specified. In addition subseq needs strand position.
		
		Arguments:
			seq 	- A list of positions for id, start and end in header
			subseq	- A list of positions for id, start, end and strand in header
		"""
		seqKeys = ['id', 'start', 'end', 'strand']
		i = 0
		for  p in seq:
			self.seqs.setPos(seqKeys[i], p)
			i += 1
		i = 0
		for p in subseq:
			self.subseqs.setPos(seqKeys[i], p)
			i += 1
		return self
		
	# Set how many letters to include before start site of union/complement
	def setUpstream(self, upstream):
		"""
		Set how many letters to include before start site of union/complent
		
		Arguments:
			upstream - Number of letters before start of sequence
		"""
		self.upstream = upstream
		return self
	
	# Option -ds
	def setDownstream(self, downstream):
		"""
		Set how many letters to include after end site of union/complement
		
		Arguments:
			downstream - Number of letters beyond end of sequence
		"""
		self.downstream = downstream
		return self
	
	# Option -g
	def setGenes(self, genes=False):
		"""
		Specify which genes to include via their ID given in header.
		
		Arguments:
			genes	- When False include all, otherwise a list of gene
					IDs to include.
		"""
		self.selectedGenes = genes
		return self
	
	# Option -u
	def setUnion(self, union=False):
		"""
		Specify which unions to include in output.
		
		Arguments:
			union	- When False include none, when True include all
						when list include all unions with order
						number in list.
		"""
		if union == []: self.getUnions = True
		else: self.getUnions = union
		return self
	
	# Option -c
	def setComplement(self, complement=False):
		"""
		Specify which complements to include.
		
		Arguments:
			complement	- When False include none, when True include all
							when list include all complements with order
							number in list.
		"""
		if complement == []: self.getComplements = True
		else: self.getComplements = complement
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
	
	# Iterate through all sequences
	def _parseGenes(self):
		o = ''
		for fname,gene in self.seqs._parse():
			try:
				metaParts = gene.id.split('|')
				geneName = metaParts[self.seqs.id]
				if self.selectedGenes == False or geneName in self.selectedGenes:
					geneStart = int(metaParts[self.seqs.start])
					geneEnd = int(metaParts[self.seqs.end])
					#self.outfile.write('>'+geneName + '\n' + self.toFASTA(str(gene.seq)) +'\n')
					exons,introns = self.subseqs.regions(geneName)
					strand = self.subseqs.positions[geneName]
					# Set sorting: When gene is on minus strand, sorting needs to be descending
					rev = False
					if strand == -1: rev = True
					if self.getUnions:
						exons.sort(reverse=rev)
						o += self._extractor(geneName, exons, 'EXON', str(gene.seq), geneStart, geneEnd, strand)
					if self.getComplements:
						introns.sort(reverse=rev)
						# If gene is on the minus strand, intron start will need to be shifted 1 position forward
						o += self._extractor(geneName, introns, 'INTRON', str(gene.seq), geneStart, geneEnd+rev, strand)
			except:
				print('Problem parsing gene with header '+gene.id+' in file '+fname)
				continue
		return o
	
	# Iterate through regions and extract them
	def _extractor(self, name, regions, type, seq, geneStart, geneEnd, strand):
		o = ''
		p = 1
		if regions != []:
			for a,b in regions:
				if (type=='INTRON' and (self.getComplements==True or p in self.getComplements)) or (type=='EXON' and (self.getUnions==True or p in self.getUnions)):
					o += '>' + name + '|' + type + '_' + str(p) + '\n' + self._toFASTA(self._getSeq(seq, geneStart, geneEnd, a, b, strand))
				p += 1
		return o
	
	# Convert coordinates to their string counterpart and return this section
	# of the sequence.
	def _getSeq(self, seq, geneStart, geneEnd, start, end, strand=1):
		# When on the minus strand the coordinates in the genome refer to basepairs
		# but the sequence given is always in the 5'->3' direction. In this case the
		# end coordinate is the starting point.
		if strand == -1:
			a = geneEnd - end
			b = geneEnd - start + 1
		else:
			a = start - geneStart
			b = end - geneStart + 1
		return seq[a-self.upstream:b+self.downstream]
	
	def __str__(self):
		self.subseqs.parseMeta()
		return self._parseGenes()
	
def main():
	# Set program description
	parser = argparse.ArgumentParser(description="""MergeSeq takes nucleotide sequences and define unions of subsequences and their complement.
	The program expects sequences in file(s) provided by -s and coordinates of subsequences in the headers of -ss. The subsequences themselves are
	not used in defining unions/complements. By default the sequence files are expected to have headers specifying ID of sequence, start position
	and end position: "ID|START|END". In addition, subsequences need to specify gene strand: "ID|START|END|STRAND". If headers have a different
	format these need to be specified with options -hs for sequences and -hss for subsequences.
	 2015 (c) Lord Jangmo Software Curp.""")
	# Set options
	parser.add_argument('-r', help='Display process report.', action='store_const', const=True, default=False)
	parser.add_argument('-f', help='Set format of input data, can be any accepted by BioPython, defaults to FASTA', default='fasta')
	parser.add_argument('-us', type=int, default=0, help='Set upstream length')
	parser.add_argument('-ds', type=int, default=0, help='Set downstream length')
	parser.add_argument('-hs', default=[1,2,3], nargs=3, type=int, help='Specify header structure in sequence file: #ID #START #END')
	parser.add_argument('-hss', default=[1,2,3,4], nargs=4, type=int, help='Specify header structure in subsequence file: #ID #START #END #STRAND')
	parser.add_argument('-g', default=False, nargs='*', help='Select gene(s) with ID')
	parser.add_argument('-u', help='Extract union of subsequences', type=int, default=False, nargs='*')
	parser.add_argument('-c', help='Extract complement to the union of subsequences', type=int, default=False, nargs='*')
	parser.add_argument('-o', type=argparse.FileType('w'), default=False, help='Direct output to file')
	parser.add_argument('-s', type=argparse.FileType('r'), nargs='+', required=True, help='File containing sequences')
	parser.add_argument('-ss', type=argparse.FileType('r'), nargs='+', required=True, help='File containing subsequences')
	try:
		# Parse command line arguments
		args = parser.parse_args()
		c = MergeSeq(args.s, args.ss).setUpstream(args.us).setDownstream(args.ds).setGenes(args.g).setUnion(args.u).setComplement(args.c).setStructure(args.hs, args.hss)
		# Direct output to file if -o is specified, otherwise stdout
		if args.o != False: args.o.write(str(c))
		else: print c,
	except KeyboardInterrupt:
		print("\nInterrupted by user.\nGoodbye!")
	except Exception as e:
		print(e)

if __name__ == '__main__':
	main()
