from argparse import ArgumentParser
import gffutils
import os
from collections import defaultdict

def get_gtf_db(gtf, in_memory=False):
    """
    create a gffutils DB from a GTF file and will use an existing gffutils
    database if it is named {gtf}.db
    """
    db_file = ":memory:" if in_memory else gtf + ".db"
    if in_memory or not os.path.exists(db_file):
        db = gffutils.create_db(gtf, dbfn=db_file)
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

if __name__ == "__main__":
    description = ("Calculate number of transcripts for each gene.")
    parser = ArgumentParser(description)
    parser.add_argument("gtf", help="GTF file to calculate")

    args = parser.parse_args()

    db = get_gtf_db(args.gtf)
    gene_dict = get_transcripts_in_genes(db)
    for gene, transcripts in gene_dict.items():
        print gene, len(transcripts)
