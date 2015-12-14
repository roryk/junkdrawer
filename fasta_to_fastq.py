"""
Convert FASTA to FASTQ file
"""

import os
import sys
import gzip
from Bio import SeqIO
from argparse import ArgumentParser

def open_fasta(in_file):
    """ open a fasta file for reading using gzip if it is gzipped
    """
    _, ext = os.path.splitext(in_file)
    if ext == ".gz":
        return gzip.open(in_file, 'rb')
    if ext in [".fasta", ".fa"]:
        return open(in_file, 'r')
    # default to just opening it
    return open(in_file, "r")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("fasta", help="path to FASTA file")
    args = parser.parse_args()

    with open_fasta(args.fasta) as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            record.letter_annotations["phred_quality"] = [60] * len(record)
            SeqIO.write(record, sys.stdout, "fastq-sanger")
