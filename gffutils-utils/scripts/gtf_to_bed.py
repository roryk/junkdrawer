import tempfile
from argparse import ArgumentParser
import gffutils
import os

def disable_infer_extent(gtf_file):
    """
    guess if we need to use the gene extent option when making a gffutils
    database by making a tiny database of 1000 lines from the original
    GTF and looking for all of the features
    """
    _, ext = os.path.splitext(gtf_file)
    tmp_out = tempfile.NamedTemporaryFile(suffix=".gtf", delete=False).name
    with open(tmp_out, "w") as out_handle:
        count = 0
        in_handle = open(gtf_file) if ext != ".gz" else gzip.open(gtf_file)
        for line in in_handle:
            if count > 1000:
                break
            out_handle.write(line)
            count += 1
        in_handle.close()
    db = gffutils.create_db(tmp_out, dbfn=":memory:",
                            disable_infer_transcripts=False,
                            disable_infer_genes=False)
    os.remove(tmp_out)
    features = [x for x in db.featuretypes()]
    if "gene" in features and "transcript" in features:
        return True
    else:
        return False

def get_gtf_db(gtf, in_memory=False):
    """
    create a gffutils DB from a GTF file and will use an existing gffutils
    database if it is named {gtf}.db
    """
    db_file = ":memory:" if in_memory else gtf + ".db"
    disable_infer = disable_infer_extent(gtf)
    if in_memory or not os.path.exists(db_file):
        db = gffutils.create_db(gtf, dbfn=db_file,
                                disable_infer_transcripts=disable_infer,
                                disable_infer_genes=disable_infer)
    if in_memory:
        return db
    else:
        return gffutils.FeatureDB(db_file)

def gtf_to_bed(gtf):
    db = get_gtf_db(gtf)
    out_file = os.path.splitext(gtf)[0] + ".bed"
    if os.path.exists(out_file):
        return out_file
    with open(out_file, "w") as out_handle:
        for feature in db.all_features():
            chrom = feature.chrom
            start = feature.start
            end = feature.end
            attributes = feature.attributes.keys()
            strand = feature.strand
            name = (feature['gene_name'][0] if 'gene_name' in attributes else
                    feature['gene_id'][0])
            line = "\t".join(map(str, [chrom, start, end, name, ".", strand]))
            out_handle.write(line + "\n")
    return out_file

if __name__ == "__main__":
    description = ("Convert a GTF file to a BED file.")
    parser = ArgumentParser(description)
    parser.add_argument("gtf", help="GTF to convert")
    args = parser.parse_args()

    gtf_to_bed(args.gtf)
