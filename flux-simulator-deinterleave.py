
import os
import sys
from Bio import SeqIO


def deinterleave(in_fastq):
    """
    removes reads from a pair of fastq files that are shorter than
    a minimum length. removes both ends of a read if one end falls
    below the threshold while maintaining the order of the reads

    """
    quality_format = "fastq-sanger"
    fq1_out = os.path.splitext(in_fastq)[0] + "_1.fq"
    fq2_out = os.path.splitext(in_fastq)[0] + "_2.fq"

    fastq_reader = SeqIO.parse(in_fastq, quality_format)

    out_files = [fq1_out, fq2_out]

    fq1_out_handle = open(out_files[0], "w")
    fq2_out_handle = open(out_files[1], "w")

    for fq1_record in fastq_reader:
        fq2_record = fastq_reader.next()
        if not fq1_record.id[:-2] == fq2_record.id[:-2]:
            continue
        fq1_out_handle.write(fq1_record.format(quality_format))
        fq2_out_handle.write(fq2_record.format(quality_format))

    return [fq1_out, fq2_out]


if __name__ == "__main__":
    deinterleave(sys.argv[1])
