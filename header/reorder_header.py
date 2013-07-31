# reorder a hg19 header into a GRCh37 header
import sys

def get_sequences_in_header(header):
    with open(header) as header_handle:
        sequences = [_sequence_from_line(x) for x in header_handle if "@SQ" in x]
    return sequences


def _sequence_from_line(line):
    return line.split("\t")[1].split(":")[1]


def _hg19_to_GRCh37(item):
    if "hap" in item:
        return None
    if "_" in item and "gl" in item:
        return item.split("_")[1].upper() + ".1"
    elif "chrM" in item:
        return "MT"
    elif "chr" in item:
        return item.replace("chr", "")


def get_hg19_converter(hg19):
    GRCh37_converted = map(_hg19_to_GRCh37, get_sequences_in_header(hg19))
    return dict(zip(get_sequences_in_header(hg19), GRCh37_converted))


def main(hg19, GRCh37):
    # hg19_sequences = get_sequences_in_header(hg19)
    # for sequence in hg19_sequences:
    #     print _hg19_to_GRCh37(sequence)
    GRCh37_sequences = get_sequences_in_header(GRCh37)
    hg19_to_GRCh37 = {}
    for item in get_sequences_in_header(hg19):
        converted_item = _hg19_to_GRCh37(item)
        if converted_item in GRCh37_sequences:
            hg19_to_GRCh37[item] = converted_item
    GRCh37_to_hg19 = dict(zip(hg19_to_GRCh37.values(), hg19_to_GRCh37.keys()))

    with open(GRCh37) as in_handle:
        for line in in_handle:
            if not "@SQ" in line:
                sys.stdout.write(line)
                continue
            line_split = line.split("\t")
            line_split[1] = "SN:" + GRCh37_to_hg19[_sequence_from_line(line)]
            sys.stdout.write("\t".join(line_split))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
