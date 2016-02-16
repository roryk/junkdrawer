from argparse import ArgumentParser
from Bio.Seq import Seq
from itertools import product

def rc_prepper(barcode):
    return str(Seq(barcode.strip()).reverse_complement())

def normal_prepper(barcode):
    return barcode.strip()

def read_barcode(fn, prepper=normal_prepper):
    with open(fn) as fh:
        return set(prepper(barcode) for barcode in fh)

if __name__ == "__main__":
    description = ("Takes in a set of barcodes and prints out their reverse "
                   "complement. If paired barcodes are given, it prints all "
                   "combinations of BC1+BC2.")
    parser = ArgumentParser()
    parser.add_argument("--rc", action='store_true', default=False,
                        help="Take reverse complement of barcodes.")
    parser.add_argument("--bc1", help="Barcode 1", required=True)
    parser.add_argument("--bc2", help="Barcode 2", default=None, required=True)

    args = parser.parse_args()

    prepper = rc_prepper if args.rc else normal_prepper

    bc1_set = read_barcode(args.bc1, prepper)
    bc2_set = read_barcode(args.bc2, prepper) if args.bc2 else None

    if bc2_set:
        for bc1, bc2 in product(bc1_set, bc2_set):
            print bc1 + bc2
    else:
        for bc1 in bc1_set:
            print bc1
