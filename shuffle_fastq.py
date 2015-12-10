# shuffles the sequences in a fastq file
import os
import random
from Bio import SeqIO
import fileinput
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--fq1", required="True")
    parser.add_argument("--fq2", required="True")
    args = parser.parse_args()
    with open(args.fq1) as in_handle:
        fq1 = [x for x in SeqIO.parse(in_handle, "fastq-sanger")]
    with open(args.fq2) as in_handle:
        fq2 = [x for x in SeqIO.parse(in_handle, "fastq-sanger")]
    order = range(len(fq1))
    random.shuffle(order)

    fq1_name = os.path.splitext(args.fq1)[0]
    fq2_name = os.path.splitext(args.fq2)[0]
    with open(fq1_name + ".shuffled.fq", "wa") as fq1_handle, open(fq2_name + ".shuffled.fq", "wa") as fq2_handle:
        for i in order:
            fq1_handle.write(fq1[i].format("fastq-sanger"))
            fq2_handle.write(fq2[i].format("fastq-sanger"))
