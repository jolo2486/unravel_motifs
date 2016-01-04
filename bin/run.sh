#!/usr/bin/env bash

# This script reproduces sequence logos

# If input files have a different name change the variables here:
# coding_sequences.fa - gene coding sequences including 15 nucleotides upstream
cs=seq/coding_sequences.fa
# usg.fa - unspliced gene sequences
usg=seq/usg.fa
# exons.fa - exon sequences
exons=seq/exons.fa


# CODING SEQUENCES - Sequence logo around the start site
# Input file $cs should be a FASTA formatted file containing coding sequences
# from all human chromosomes including 15 characters upstream from the
# start codon. The cutseq command extracts the first 33 nucleotides in each
# sequence in $cs, thus the start codon (most common ATG) is aligned at
# position 16 in the extracted sequences.
echo "CODING SEQUENCES"
echo "Extracting first 33 characters from coding sequences, including 15 nucleotides upstream..."
./cutseq.py -l 33 $cs > output/coding_sequences_1_33.fa
echo "...done!"
# Create the logo with the createlogo command, y-scale is set to 0.3 bits
# to make the motifs around the start site more visible.
echo "Creating logo..."
./createlogo.py -y 0.3 -o output/logo_coding_sequences.png output/coding_sequences_1_33.fa
echo "...done!"

# INTRON SEQUENCES - Sequence logos around start and end of first intron in each gene
# Input file $usg should contain unspliced genes from all human chromosomes
# in FASTA format where headers should have the format
# "GENE ID|GENE START POSITION|GENE END POSITION" (without having to modify
# any command line options). Input file $exons should contain exon sequences
# from all human genes with headers specified as "GENE ID|EXON START|EXON END|STRAND".
# The command mergeseq takes $usg and $exons and with option "-c 1" the first
# intron in each gene is extracted by first taking the complement of all exon
# sequences specified by coordinates and strand in $exons.
echo "INTRON SEQUENCES"
echo "Extracting first intron, 15 nucleotides upstream/downstream..."
./mergeseq.py -us 15 -ds 15 -c 1 -s $usg -ss $exons -o output/usg_intr_1.fa
echo "..done!"
echo "Extracting first 30 nucleotides from intron..."
./cutseq.py -l 30 -o output/usg_intr_1_us_15.fa output/usg_intr_1.fa
echo "...done!"
echo "Extracting last 30 nucleotides from intron..."
./cutseq.py -l -30 -o output/usg_intr_1_ds_15.fa output/usg_intr_1.fa
echo "...done!"
# Create upstream logo
echo "Creating intron upstream logo..."
./createlogo.py -y 2 -o output/logo_usg_intr_1_us_15.png output/usg_intr_1_us_15.fa
echo "...done!"
# Create downstream logo
echo "Creating intron downstream logo..."
./createlogo.py -o output/logo_usg_intr_1_ds_15.png output/usg_intr_1_ds_15.fa
echo "...done!"
