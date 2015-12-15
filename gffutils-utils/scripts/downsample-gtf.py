"""
Downsample a GTF file, dropping random transcripts from genes.
"""

from argparse import ArgumentParser
import gffutils
import os
from random import sample
from collections import defaultdict

def get_gtf_db(gtf, in_memory=False):
    """
    create a gffutils DB from a GTF file and will use an existing gffutils
    database if it is named {gtf}.db
    """
    db_file = ":memory:" if in_memory else gtf + ".db"
    if in_memory or not os.path.exists(db_file):
        db = gffutils.create_db(gtf, dbfn=db_file, disable_infer_transcripts=True)
    if in_memory:
        return db
    else:
        return gffutils.FeatureDB(db_file)

def get_transcripts_in_genes(db):
    """
    return a dictionary where the keys are gene_ids and the values are the
    set of transcript_ids assigned to that gene_id
    """
    genes = defaultdict(set)
    for feature in db.all_features():
        gene_id = feature.attributes.get('gene_id', [None])[0]
        transcript_id = feature.attributes.get('transcript_id', [None])[0]
        if gene_id and transcript_id:
            genes[gene_id].update([transcript_id])
    return genes


def get_longest_transcripts(db):
    """
    returns two lists, the first is the longest transcripts and the
    second are the rest of the transcripts
    """
    longest = []
    transcripts = []
    for gene in db.features_of_type('gene'):
        max_length = 0
        max_transcript = ''
        for transcript in db.children(gene, level=1):
            transcripts.append(transcript)
            exon_lengths = []
            for exon in db.children(transcript, level=1, featuretype='exon'):
                exon_lengths.append(abs(exon.start - exon.end))
                if sum(exon_lengths) > max_length:
                    max_transcript = transcript
                    max_length = sum(exon_lengths)
        longest.append(transcript)
    transcripts = [x for x in transcripts if x not in longest]
    return longest, transcripts

if __name__ == "__main__":
    description = "Downsample a GTF file, dropping random transcripts from genes."
    parser = ArgumentParser(description=description,
                            usage=__doc__.strip())
    parser.add_argument("--percentage", default=0.5,
                        help="Percentage of transcripts to keep.")
    parser.add_argument("--keep-longest", action="store_true",
                        default=False,
                        help="Never drop the longest transcript for each gene")
    parser.add_argument("gtf", help="GTF file to calculate")

    args = parser.parse_args()

    db = get_gtf_db(args.gtf)
    longest, transcripts = get_longest_transcripts(db)
    if args.keep_longest:
        for transcript in longest:
            print transcript
    else:
        transcripts = transcripts + longest
    keep = sample(transcripts, int(len(transcripts) * args.percentage))
    for k in keep:
        print k
