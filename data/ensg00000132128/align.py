#!/usr/bin/env python2.7
from Bio import SeqIO
import math,time
startTime = time.time()

# Pad string with a character to desired width
def pad(seq, start, geneLen, pad=' '):
	endAdd = geneLen - start - len(seq)
	return pad * start + seq + pad * endAdd + '\n'

# Convert string to FASTA format
def toFASTA(seq, width=60):
	o = ''
	a = 0
	b = 0
	for x in range(1, int(math.floor(len(seq)/60))):
		b = 60 * x
		o += seq[a:b] + '\n'
		a = b
	return o + seq[a+60:]
	
# Files for writing
# Plain alignment 5'|-------|3'
f = open('alignment.txt', 'w')
# Fasta formatted
fasta = open('alignment.fa', 'w')
# Load gene
usg = list(SeqIO.parse('usg.fa', 'fasta'))[0]
gId = str(usg.id)
gseq = str(usg.seq)
gSeqLen = len(gseq)
# Output unspliced gene
f.write(gId[0:15] + ' ' + gseq + '\n')
fasta.write('>' + gId[0:15] + '\n' + toFASTA(gseq) +'\n')
# Now iterate through transcripts
positions = []
for ts in SeqIO.parse('ust.fa', 'fasta'):
	# Transcripts
	tseq = str(ts.seq)
	tId = str(ts.id)[0:15]
	pos = gseq.find(tseq)
	if pos == -1: pos = 0
	fasta.write('>---' + tId[3:] + '\n' + toFASTA(pad(tseq, pos, gSeqLen, '-')))
	f.write('   ' + tId[3:] + ' ' + pad(tseq, pos, gSeqLen, '.'))
	# 5'-UTR
	for utr in SeqIO.parse('5utr.fa', 'fasta'):
		if str(utr.id) == tId:
			u5Seq = str(utr.seq)
			pos = gseq.find(u5Seq)
			if pos == -1: utrs = '\n'
			else: utrs = pad(u5Seq, pos, gSeqLen, '-')
			fasta.write('>'+'-'*11+'UTR5\n'+toFASTA(utrs))
			f.write(' '*11 + 'UTR5 ' + utrs)
	# Exons
	for ex in SeqIO.parse('exons.txt', 'fasta'):
		exTId = str(ex.id).split('|')[0]
		exId = str(ex.id).split('|')[-1]
		if tId == exTId:
			exSeq = str(ex.seq)
			pos = gseq.find(exSeq)
			positions += [(pos, pos+len(exSeq))]
			if pos == -1: pos = 0
			fasta.write('>---'+ exId[3:] + '\n' + toFASTA(pad(exSeq, pos, gSeqLen, '-')))
			f.write('   '+ exId[3:] + ' ' + pad(exSeq, pos, gSeqLen, '.'))
	# 3'-UTR
	for utr in SeqIO.parse('3utr.fa', 'fasta'):
		if str(utr.id) == tId:
			u3Seq = str(utr.seq)
			pos = gseq.find(u3Seq)
			if pos == -1: utrs = '\n'
			else: utrs = pad(u3Seq, pos, gSeqLen, '-')
			fasta.write('>'+'-'*11+'UTR3\n'+toFASTA(utrs))
			f.write(' '*11 + 'UTR3 ' + utrs)
f.close()
positions.sort()
a = -1
b = False
# Identify the union of exons
exons = []
for x in positions:
	if a == -1 or x[0] < a: a = x[0]
	if x[0] > b:
		exons += [(a,b)]
		a = x[0]
		b = x[1]
	if b == False or x[1] > b: b = x[1]
if exons[-1][1] < a: exons += [(a,b)]
# The introns are now the complement to the exons
i = 1
p = 0
for x in exons[1:]:
	a = exons[p][1]
	b = x[0] - a
	seq = gseq[a:a+b]
	fasta.write('>INTRON_'+str(i)+'\n'+toFASTA(pad(seq, a, gSeqLen, '-')))
	i += 1
	p += 1
print exons
fasta.close()
print time.time() - startTime
