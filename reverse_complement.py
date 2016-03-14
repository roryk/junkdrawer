import sys
from argparse import ArgumentParser, FileType
from Bio.Seq import Seq

def reverse_complement(seq):
    return str(Seq(seq.strip()).reverse_complement())

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("infile", nargs="?", type=FileType('r'), default=sys.stdin)
    args = parser.parse_args()

    for seq in args.infile:
        print reverse_complement(seq)
