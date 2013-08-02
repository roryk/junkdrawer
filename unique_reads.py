# filter a BAM file to extract only the uniquely mapped reads
import pysam
import sys
from argparse import ArgumentParser


def main(in_bam, out_bam, max_multi=1):
    samfile = pysam.Samfile(in_bam, "rb")
    outfile = pysam.Samfile(out_bam, "wb", template=samfile)
    for read in samfile.fetch():
        for tag in read.tags:
            if tag[0] == 'NH':
                if int(tag[1]) <= max_multi:
                    outfile.write(read)
    samfile.close()
    outfile.close()

if __name__ == "__main__":
    parser = ArgumentParser(description="Filter reads mapping more than "
                            "max_multi times out of a BAM file.")
    parser.add_argument("--max-multi", default=1, type=int,
                        help="Maximum number of multi-hits to allow per read.")
    parser.add_argument('in_bam', help="Input BAM file name.")
    parser.add_argument('out_bam', help="Output BAM file name.")
    args = parser.parse_args()
    main(args.in_bam, args.out_bam, args.max_multi)
