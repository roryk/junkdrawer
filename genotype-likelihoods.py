from __future__ import print_function
import sys
import cyvcf
from argparse import ArgumentParser, FileType
import toolz as tz

description = ("Create a table of probability of a non reference call for each "
               "genotype for each sample. This is PL[0]. -1 is output for samples "
               "with a missing PL call at a position.")
parser = ArgumentParser(description=description)
parser.add_argument("vcf", type=FileType('r'),
                    help="VCF file to convert, use '-' to read from stdin")
args = parser.parse_args()

vcf_reader = cyvcf.Reader(args.vcf)

samples = vcf_reader.samples[1:5]

header = "\t".join([str(x) for x in ["CHROM", "POS", "ID", "REF", "ALT"] + samples])

print(header, file=sys.stdout)
for record in vcf_reader:
    line = [record.CHROM, record.POS, record.ID, record.REF, record.alleles[1]]
    pls = [x.data.get("PL", None) for x in record.samples[1:5]]
    pls = [x[0] if x else "-1" for x in pls]
    print("\t".join([str(x) for x in line + pls]), file=sys.stdout)
