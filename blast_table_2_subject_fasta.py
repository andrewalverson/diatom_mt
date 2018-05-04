#! /usr/bin/env python3

# load required modules
import sys
import os
import re
import argparse
from Bio import SeqIO
from collections import defaultdict


def get_args():
	# create an ArgumentParser object ('parser')
	parser = argparse.ArgumentParser(description="This script reads a BLAST table and calls eutils to extract subject matches from NCBI")

	# add required/positional arguments
	parser.add_argument("blast",  help="BLAST output: -outfmt '6 std qcovs'")
	parser.add_argument("-e", "--emax",  help="maximum e-value", default=1e-12, type=float)
	parser.add_argument("-q", "--qcovs", help="maximum e-value", default=50, type=int)

	# parse the arguments
	return parser.parse_args()


def parse_blast():
	# open and read BLAST output
	blast = open(args.blast, 'r')

	# key = accession number, value = list of sstart + send hsp coords; will take min and max coords to extract the "full" match when there are multiple HSPs
	blast_coords = defaultdict(dict) 
	# key = accession number, value = strand
	strands = defaultdict(dict) 

	for line in blast:
		# parse BLAST line
		(query, sbjct, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qcovs) = line.rstrip().split('\t')

		if((int(qcovs) < args.qcovs) or (float(evalue) > args.emax)): 
			continue

		# calculate and store the strand
		if(int(sstart) < int(send)):
			strands[sbjct] = 1  # '+' strand
		else:
			strands[sbjct] = 2  # '-' strand

		# add start and end coords to a list of coords for this subject sequence
		if sbjct in blast_coords:
			# print("in:", sbjct, sstart, send, evalue, qcovs)
			blast_coords[sbjct].extend([int(sstart), int(send)])
		else:
			# print("out:", sbjct, sstart, send, evalue, qcovs)
			blast_coords[sbjct] = []
			blast_coords[sbjct].extend([int(sstart), int(send)])


	blast.close()

	return(blast_coords, strands)


def fetch_seqs(blast_coords, strands):
	for sbjct in blast_coords.keys():
		# assemble esearch + efetch command
		esearch = "esearch" + " -db nuccore" + " -query " + sbjct + " |" + " efetch -format fasta" + " -seq_start " + str(min(blast_coords[sbjct])) + " -seq_stop " + str(max(blast_coords[sbjct])) + " -strand " + str(strands[sbjct])

		# execute esearch command
		# print(esearch)
		os.system(esearch)


def main():
	# parse the BLAST output
	blast_coords, strands = parse_blast()

	# fetch and print the subject sequences
	fetch_seqs(blast_coords, strands)

# get the command-line arguments
args = get_args()

# execute the program by calling main
if __name__ == '__main__':
	main()
