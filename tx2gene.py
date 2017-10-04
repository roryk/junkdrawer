from __future__ import print_function
from argparse import ArgumentParser

def strip_non_alphanumeric(string):
    return string.replace('"', '').replace(';', '')

def tx2featuredict(gtf, feature="gene_id"):
    """
    produce a tx2gene dictionary from a GTF file. assumes the GTF has
    feature and transcript_id set at a minimum
    """
    d = {}
    with open(gtf) as in_handle:
        for line in in_handle:
            if not (feature in line and "transcript_id" in line):
                continue
            featureid = line.split("gene_id")[1].split(" ")[1]
            featureid = strip_non_alphanumeric(featureid)
            txid = line.split("transcript_id")[1].split(" ")[1]
            txid = strip_non_alphanumeric(txid)
            d[txid] = featureid
    return d

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("gtf", help="GTF file")
    parser.add_argument("--feature", default="gene_id",
                        help="which GTF feature to map transcripts to")
    args = parser.parse_args()
    for k, v in tx2featuredict(args.gtf, args.feature).items():
	    print(",".join([k, v]))
