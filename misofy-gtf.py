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
import logging
from argparse import ArgumentParser

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s %(name)s %(levelname)s %(message)s',
                    datefmt='%m-%d-%y %H:%M:%S')
logger = logging.getLogger("misofy-gtf")

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

def make_miso_annotation(ucsc, outdir, genomename):
    cmd = ("python gff_make_annotation.py --genome-label {genomename} {outdir} {outdir}")
    subprocess.check_call(cmd.format(**locals()), shell=True)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("gtf", help="GTF file to convert")
    parser.add_argument("outdir", help="output directory")
    parser.add_argument("genomename", help="genome name")
    args = parser.parse_args()

    gtf = args.gtf
    logger.info("Starting conversion of %s to MISO event annotation." % gtf)
    outdir = args.outdir
    genomename = args.genomename
    safe_makedir(outdir)
    logger.info("Converting %s to genepred format." % gtf)
    genepred = gtf_to_genepred(gtf, outdir)
    logger.info("Converting %s to UCSC table format." % genepred)
    ucsc = genepred_to_UCSC_table(genepred, outdir)
    logger.info("Making MISO event annotation.")
    make_miso_annotation(ucsc, outdir, genomename)
    logger.info("Finished.")
