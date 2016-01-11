#!/usr/bin/env bash

# This script reproduces sequence logos

# INSTRUCTIONS FOR REPRODUCING THE RESULTS

# FILES
# The gene data has to be downloaded manually due to file size restrictions on github.
# The queries for the files can be found in data/human_genome_biomart_sequence_links.txt.
# the files will have to be renamed to coding_sequences.fa, usg.fa, exons.fa (in the same 
# order as they appear in the above mentioned query file). They should be approx 1.4Gb in total.

# Note: Biomart seems to have some issues regarding URL queries like this. This has unfortunately
# come to our attention very late in the project so to fix it we recommend that, in events of the
# script not working, you type in the queries manually at Biomart as they are described in
# human_genome_biomart_sequence_links.txt.

# The files are to be put in a directory called seq directly under the main project. If other
# location is used, please change the variables cs, usg, exons, further down.

# OUTPUT
# The results will be generated and put in results/output. This will include some temporary
# sequence fasta files that our program uses. These could easily have been deleted by the
# script but are left for inspection purposes.

# IN SHORT
# - Download the whole project from github: unravel_motifs/
# - Download the sequence files using queries from data/human_genome_biomart_sequence_links.txt.
# - Create unravel_motifs/seq, rename and put the sequence files here.
# - Run this script: unravel_motifs/bin/run.sh
# - Results will be created in unravel_motifs/results/output

# If input files have a different name change the variables here:
# coding_sequences.fa - gene coding sequences including 15 nucleotides upstream
cs=../seq/coding_sequences.fa
# usg.fa - unspliced gene sequences
usg=../seq/usg.fa
# exons.fa - exon sequences
exons=../seq/exons.fa
# output directory
output=../results/output

# CODING SEQUENCES - Sequence logo around the start site
# Input file $cs should be a FASTA formatted file containing coding sequences
# from all human chromosomes including 15 characters upstream from the
# start codon. The cutseq command extracts the first 33 nucleotides in each
# sequence in $cs, thus the start codon (most common ATG) is aligned at
# position 16 in the extracted sequences.
mkdir ../results/output
echo "CODING SEQUENCES"
echo "Extracting first 33 characters from coding sequences, including 15 nucleotides upstream..."
./cutseq.py -l 33 $cs > $output/coding_sequences_1_33.fa
echo "...done!"
# Create the logo with the createlogo command, y-scale is set to 0.3 bits
# to make the motifs around the start site more visible.
echo "Creating logo..."
./createlogo.py -y 0.3 -o $output/logo_coding_sequences.png $output/coding_sequences_1_33.fa
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
./mergeseq.py -us 15 -ds 15 -c 1 -s $usg -ss $exons -o $output/usg_intr_1.fa
echo "..done!"
echo "Extracting first 30 nucleotides from intron..."
./cutseq.py -l 30 -o $output/usg_intr_1_us_15.fa $output/usg_intr_1.fa
echo "...done!"
echo "Extracting last 30 nucleotides from intron..."
./cutseq.py -l -30 -o $output/usg_intr_1_ds_15.fa $output/usg_intr_1.fa
echo "...done!"
# Create upstream logo
echo "Creating intron upstream logo..."
./createlogo.py -y 2 -o $output/logo_usg_intr_1_us_15.png $output/usg_intr_1_us_15.fa
echo "...done!"
# Create downstream logo
echo "Creating intron downstream logo..."
./createlogo.py -o $output/logo_usg_intr_1_ds_15.png $output/usg_intr_1_ds_15.fa
echo "...done!"
