"""
write MISO summary tables out to tidy format
"""

from argparse import ArgumentParser
import pandas as pd
import os
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(name)s %(levelname)s %(message)s',
                    datefmt='%m-%d-%y %H:%M:%S')
logger = logging.getLogger("miso-compare")

def is_misosummary(filename):
    return filename.endswith(".miso_summary")

def get_summarytype(filename):
    stem = os.path.basename(os.path.splitext(filename)[0])
    return stem.split("-")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("summarydir")
    args = parser.parse_args()
    outfile = os.path.join(args.summarydir, "combined-miso.csv")

    misofiles = [x for x in os.listdir(args.summarydir) if is_misosummary(x)]
    frame = pd.DataFrame()
    reslist = []
    for misofile in misofiles:
        logger.info("Parsing %s." % misofile)
        samplename, eventtype = get_summarytype(misofile)
        misopath = os.path.join(args.summarydir, misofile)
        df = pd.read_table(misopath, sep="\t", header=0)
        df['samplename'] = samplename
        df['eventtype'] = eventtype
        reslist.append(df)
    frame = pd.concat(reslist)
    logger.info("Writing tidy MISO summaries to %s." % outfile)
    frame.to_csv(outfile, index=False)
