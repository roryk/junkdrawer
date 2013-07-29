# filter a BAM file to extract only the uniquely mapped reads
import pysam
import sys


def main(in_bam, out_bam):
    samfile = pysam.Samfile(in_bam, "rb")
    outfile = pysam.Samfile(out_bam, "wb", template=samfile)
    for read in samfile.fetch():
        if ('NH', 1) in read.tags:
            outfile.write(read)
    samfile.close()
    outfile.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
