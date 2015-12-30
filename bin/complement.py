#!/usr/bin/env python2.7
from Bio import SeqIO
import math,time

a = time.time()
baseDir = '../data/ensg00000132128/'

class Sequences:
	def __init__(self, fileResources, format='fasta'):
		if isinstance(fileResources, list) == False: self.fileResources = [fileResources]
		else: self.fileResources = fileResources
		self.format = format
		self.meta = {}
		self.positions = {}
	
	def setPos(self, name, position):
		self.__dict__[name] = position - 1
		return self
	
	def parseMeta(self):
		for fr in self.fileResources:
			for seq in SeqIO.parse(fr, self.format):
				parts = seq.id.split('|')
				try:
					coord = [(int(parts[self.start]),int(parts[self.end]))]
					id = parts[self.id]
					if id in self.meta: self.meta[id] += coord
					else: self.meta[id] = coord
				except:
					1
		return self

	def regions(self, gene, relative=True):
		coordinates = self.meta[gene]
		coordinates.sort()
		exonCoord = []
		intronCoord = []
		coordinates.sort()
		start, end = coordinates[0]
		for coord in coordinates[1:]:
			if coord[0] > end:
				exonCoord += [(start,end)]
				intronCoord += [(end+1,coord[0])]
				start, end = coord
			if coord[0] < start: start = coord[0]
			if coord[1] > end: end = coord[1]
		exonCoord += [(start,end)]
		return exonCoord,intronCoord

class Complement:
	def __init__(self, geneResources, exonResources, format='fasta'):
		if isinstance(geneResources, list) == False: self.geneResources = [geneResources]
		else: self.geneResources = geneResources
		self.seqs = Sequences(exonResources).setPos('id', 1).setPos('start', 2).setPos('end', 3).parseMeta()
		self.format = format
		self.outfile = open('coord_alignment.fa', 'w')
		self._parseGenes()
	
	# Pad string with a character to desired width
	def pad(self, seq, start, geneLen, pad='-'):
		endAdd = geneLen - start - len(seq)
		return pad * start + seq + pad * endAdd + '\n'

	# Convert string to FASTA format
	def toFASTA(self, seq, width=60):
		o = ''
		a = 0
		b = 0
		for x in range(1, int(math.floor(len(seq)/60))):
			b = 60 * x
			o += seq[a:b] + '\n'
			a = b
		return o + seq[a+60:]
	
	def _parseGenes(self):
		for res in self.geneResources:
			for gene in SeqIO.parse(res, self.format):
				metaParts = gene.id.split('|')
				geneName = metaParts[0]
				geneStart = int(metaParts[-2])
				geneEnd = int(metaParts[-1])
				self.outfile.write('>'+geneName + '\n' + self.toFASTA(str(gene.seq)) +'\n')
				exons, introns = self.seqs.regions(geneName)
				introns.sort(reverse=True)
				exons.sort(reverse=True)
				p = 1
				for x in exons:
					start = geneEnd - x[1]
					end = geneEnd - x[0] + 1
					self.outfile.write('>EXON_'+str(p)+'\n'+self.toFASTA(self.pad(str(gene.seq)[start:end], start, geneEnd-geneStart+1)))
					p += 1
				p = 1
				for x in introns:
					# Intron is shifted 1 step rightwards when on -1 strand
					start = geneEnd - x[1] + 1
					end = geneEnd - x[0] + 1
					self.outfile.write('>INTRON_'+str(p)+'\n'+self.toFASTA(self.pad(str(gene.seq)[start:end], start, geneEnd-geneStart+1)))
					p += 1
	
Complement(open(baseDir+'usg.fa'), open(baseDir+'exons2.fa'))
#b = Sequences(open(baseDir + 'exons3.txt', 'r')).setPos('id', 1).setPos('start', 5).setPos('end', 6).setPos('strand', 7)
#b.parseMeta()
#print b.exons('ENSG00000132128')
print time.time()-a
