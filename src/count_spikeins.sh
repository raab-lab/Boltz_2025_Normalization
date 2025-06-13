#!/bin/bash
##########################################
# Written by Dr. Bryan Venters, EpiCypher Inc. 
# Updated 29 OCT 2021
## Updated 14 AUG 2024 By Peyton Kuhlers
## Raab Lab UNC
##########################################

# Barcode identities
# Unmodified (A & B)
# TTCGCGCGTAACGACGTACCGT 
# CGCGATACGACCGCGTTACGCG 

# H3K4me1 (A & B)
# CGACGTTAACGCGTTTCGTACG 
# CGCGACTATCGCGCGTAACGCG 

# H3K4me2 (A & B)
# CCGTACGTCGTGTCGAACGACG 
# CGATACGCGTTGGTACGCGTAA 

# H3K4me3 (A & B)
# TAGTTCGCGACACCGTTCGTCG 
# TCGACGCGTAAACGGTACGTCG

# H3K9me1 (A & B)
# TTATCGCGTCGCGACGGACGTA
# CGATCGTACGATAGCGTACCGA

# H3K9me2 (A & B)
# CGCATATCGCGTCGTACGACCG
# ACGTTCGACCGCGGTCGTACGA

# H3K9me3 (A & B)
# ACGATTCGACGATCGTCGACGA
# CGATAGTCGCGTCGCACGATCG

# H3K27me1 (A & B)
# CGCCGATTACGTGTCGCGCGTA
# ATCGTACCGCGCGTATCGGTCG

# H3K27me2 (A & B)
# CGTTCGAACGTTCGTCGACGAT
# TCGCGATTACGATGTCGCGCGA

# H3K27me3 (A & B)
# ACGCGAATCGTCGACGCGTATA
# CGCGATATCACTCGACGCGATA

# H3K36me1 (A & B)
# CGCGAAATTCGTATACGCGTCG
# CGCGATCGGTATCGGTACGCGC

# H3K36me2 (A & B)
# GTGATATCGCGTTAACGTCGCG
# TATCGCGCGAAACGACCGTTCG

# H3K36me3 (A & B)
# CCGCGCGTAATGCGCGACGTTA
# CCGCGATACGACTCGTTCGTCG

# H4K20me1 (A & B)
# GTCGCGAACTATCGTCGATTCG
# CCGCGCGTATAGTCCGAGCGTA

# H4K20me2 (A & B)
# CGATACGCCGATCGATCGTCGG
# CCGCGCGATAAGACGCGTAACG

# H4K20me3 (A & B)
# CGATTCGACGGTCGCGACCGTA
# TTTCGACGCGTCGATTCGGCGA

# template loop begin 
## Expects Raab lab style sample sheet with Read 1 in the 8th column and Read 2 in the 9th
# can only run one of the R1/R2 at a time, which is why R2 is commented out
# the 'BARCODES' is a text file that just lists the above barcode sequences (just the sequences not the names, only one colum)

SAMPLESHEET='/proj/jraablab/users/jlboltz/Normalization/norm_ss.csv'
BARCODES='/proj/jraablab/users/pkuhlers/seq_resources/spikein_barcodes.txt'

awk -F ',' 'NR>1{print $8}' $SAMPLESHEET | \ 
	parallel -j12 --tag -a $BARCODES -a - zgrep -c > /proj/jraablab/users/jlboltz/Normalization/R1_spikein_counts.txt

# awk -F ',' 'NR>1{print $9}' $SAMPLESHEET  | \
	parallel -j12 --tag -a $BARCODES -a - zgrep -c > /proj/jraablab/users/jlboltz/Normalization/R2_spikein_counts.txt

# template loop end ##
