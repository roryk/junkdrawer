# clip reads in a BAM file to a specified length
import pysam
import sys


def main(in_bam, target_length, out_bam):
    samfile = pysam.Samfile(in_bam, "rb")
    outfile = pysam.Samfile(out_bam, "wb", template=samfile)
    for read in samfile.fetch():
        if len(read.seq) < target_length:
            continue
        elif len(read.seq) > target_length:
            new_qual = read.qual[0:target_length]
            read.seq = read.seq[0:target_length]
            read.qual = new_qual
        outfile.write(read)
    samfile.close()
    outfile.close()


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
