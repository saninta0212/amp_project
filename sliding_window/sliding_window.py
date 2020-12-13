#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO

import os
import sys

# Ask for input of
# Important: 1st argument is path to FASTA, second argument is output directory, third argument is sliding window
args = sys.argv[1:]

if len(sys.argv[1:]) != 3:
    script_name = sys.argv[0]
    print(f"Usage: {script_name} <FASTA> <output directory> <size of sliding window>", file=sys.stderr)
    sys.exit(1)

INPUT_FASTA = args[0]
OUTDIR = args[1]
MER_SIZE = args[2]

if os.path.isdir(OUTDIR):
    print("Directory already exists. Proceeding with cleaving...", file=sys.stderr)
else:
    print("Directory does NOT exist. Create a new directory and try again...", file=sys.stderr)
    sys.exit(1)

cleaved_seq = open(OUTDIR + '/sliding_window_cleaved.' + MER_SIZE + 'mer.faa', "w+")

# Read FASTA input
fasta_input = open(INPUT_FASTA)
for record in SeqIO.parse(fasta_input, "fasta"):
    sequence = str(record.seq)
    n = len(sequence)

    # Cleaving using string slicing
    if n > int(MER_SIZE):
        initial = 0
        end = int(MER_SIZE)
        while end <= (n):
            cleaved_piece = sequence[initial:end]
            cleaved_seq.write(">" + str(record.id) + "_cleave_" + str(initial + 1) + "-" + str(end) + "\n" +
                              str(cleaved_piece) + "\n")
            initial = initial + 1
            end = end + 1
    else:
        cleaved_seq.write(">" + str(record.id) + "_full" + "\n" + str(sequence) + "\n")
        print(str(record.id) + ": sequence is shorter than the mer_size!", file=sys.stderr)

cleaved_seq.close()
