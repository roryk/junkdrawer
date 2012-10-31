import vcf
from collections import OrderedDict
from functools import reduce
import csv
import sys

FIELDS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]


def _make_fields_dict(record):
    return dict([(field, getattr(record, field)) for field in FIELDS])


def _record_to_dict(record):
    return OrderedDict(_make_fields_dict(record).items() +
                       record.INFO.items())


def vcf_to_tsv(in_file, out_file=None):
    if out_file is None:
        out_file = ".".join(in_file.split(".")[:-1] + ["csv"])

        writer = None

    with open(in_file, 'rb') as in_handle, open(out_file, 'w') as out_handle:
        vcf_reader = vcf.Reader(in_handle)
        header = vcf_reader.infos.keys() + FIELDS
        writer = csv.DictWriter(out_handle, fieldnames=header,
                                delimiter="\t", quotechar="/",
                                quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        for record in vcf_reader:
            d = _record_to_dict(record)
            d = {k: d.get(k, None) for k in header}
            writer.writerow(d)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        vcf_to_tsv(sys.argv[1])
    if len(sys.argv) == 3:
        vcf_to_tsv(sys.argv[1], sys.argv[2])
    else:
        print "usage: python vcf2tsv.py in_vcf [out_csv]"
