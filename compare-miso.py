"""
this script is to do all pairwise comparisons of each class of events
called by MISO. It assumes MISO has already been run on each sample
and there is a directory structure of:

miso_dir/control-RI
miso_dir/knockdown-RI

where before the - is the samplename and after the - is the event type.
It then calculates all pairwise comparisons of samples for each event type.
"""
import fnmatch
import os
import logging
import subprocess
from argparse import ArgumentParser
from itertools import combinations
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(name)s %(levelname)s %(message)s',
                    datefmt='%m-%d-%y %H:%M:%S')
logger = logging.getLogger("miso-compare")

EVENTTYPES = ["A3SS", "A5SS", "MXE", "RI", "SE"]

def is_misodir(dirname):
    for etype in EVENTTYPES:
        if dirname.endswith("-" + etype):
            return True
    return False

def drop_after_last(string, dropafter):
    """drop everything after the last match of dropafter"""
    tokens = string.split(dropafter)[:-1]
    return(dropafter.join(tokens))

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--outdir", default="comparisons")
    parser.add_argument("misodir", help="Toplevel MISO results directory.")
    args = parser.parse_args()
    misodirs = [x for x in os.listdir(args.misodir) if is_misodir(x)]
    samples = set([drop_after_last(x, "-") for x in misodirs])
    typestorun = set([x.split("-")[-1] for x in misodirs])
    pairs = list(combinations(samples, 2))
    miso_cmd = "compare_miso --compare-samples {s1} {s2} {args.outdir}"
    for etype in typestorun:
        for pair in pairs:
            s1 = pair[0] + "-" + etype
            s2 = pair[1] + "-" + etype
            cmd = miso_cmd.format(**locals())
            logger.info("Comparing %s and %s." %(s1, s2))
            subprocess.check_call(cmd, shell=True)
