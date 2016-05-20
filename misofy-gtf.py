"""
this script takes a GTF and converts it into event level files for use
with MISO. This does no sanitizing of the GTF or sanity checking or anything,
it is totally stupid, like me. If you run into problems using this post an
issue at https://github.com/roryk/junkdrawer and I'll try to help straighten
it out.

This assumes you have the following installed and in your path:
gtfToGenePred:
conda install -c bioconda ucsc-gtftogenepred

gff_make_annotation.py:
wget https://raw.githubusercontent.com/yarden/rnaseqlib/clip/rnaseqlib/gff/gff_make_annotation.py)

and

pip install misopy
"""

import os
import subprocess
from argparse import ArgumentParser

def swapdir(filepath, dirname):
    return os.path.join(dirname, os.path.basename(filepath))

def safe_makedir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def gtf_to_genepred(gtf, outdir):
    out_file = os.path.splitext(swapdir(gtf, outdir))[0] + ".genePred"
    if os.path.exists(out_file):
        return out_file
    cmd = ("gtfToGenePred -allErrors -ignoreGroupsWithoutExons -genePredExt {gtf} "
           "{out_file}")
    subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def genepred_to_UCSC_table(genepred, outdir):
    header = ["#bin", "name", "chrom", "strand",
              "txStart", "txEnd", "cdsStart", "cdsEnd",
              "exonCount", "exonStarts", "exonEnds", "score",
              "name2", "cdsStartStat", "cdsEndStat",
              "exonFrames"]
    out_file = os.path.join(outdir, "ensGene.txt")
    if os.path.exists(out_file):
        return out_file
    with open(genepred) as in_handle, open(out_file, "w") as out_handle:
        counter = -1
        current_item = None
        out_handle.write("\t".join(header) + "\n")
        for l in in_handle:
            item = l.split("\t")[0]
            if current_item != item:
                current_item = item
                counter = counter + 1
            out_handle.write("\t".join([str(counter), l]))
    return out_file

def make_miso_annotation(ucsc, outdir):
    cmd = ("python gff_make_annotation.py {outdir} {outdir}")
    subprocess.check_call(cmd.format(**locals()), shell=True)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("gtf", help="GTF file to convert")
    parser.add_argument("outdir", help="output directory")
    args = parser.parse_args()

    gtf = args.gtf
    outdir = args.outdir
    safe_makedir(outdir)
    genepred = gtf_to_genepred(gtf, outdir)
    ucsc = genepred_to_UCSC_table(genepred, outdir)
    make_miso_annotation(ucsc, outdir)
