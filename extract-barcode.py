"""
code to extract a single cell from a set of alignments or reads marked via Valentine's umis 
repository:
https://github.com/vals/umis
"""
import regex as re
import sys
from argparse import ArgumentParser
from pysam import AlignmentFile

def extract_barcode(sam, barcode):

    parser_re = re.compile('.*:CELL_(?P<CB>.*):UMI_(?P<MB>.*)')
    sam_file = AlignmentFile(sam, mode='r')
    filter_file = AlignmentFile("-", mode='wh', template=sam_file)
    track = sam_file.fetch(until_eof=True)
    for i, aln in enumerate(track):
        if aln.is_unmapped:
            continue
        match = parser_re.match(aln.qname)
        CB = match.group('CB')
        if CB == barcode:
            filter_file.write(aln)

def stream_fastq(file_handler):
    ''' Generator which gives all four lines if a fastq read as one string
    '''
    next_element = ''
    for i, line in enumerate(file_handler):
        next_element += line
        if i % 4 == 3:
            yield next_element
            next_element = ''

def extract_barcode_fastq(fastq, barcode):
    parser_re = re.compile('.*:CELL_(?P<CB>.*):UMI_(?P<MB>.*)')
    fastq_file = stream_fastq(open(fastq))
    for read in fastq_file:
        match = parser_re.match(read)
        CB = match.group('CB')
        if CB == barcode:
	    sys.stdout.write(read)

if __name__ == "__main__":
    parser = ArgumentParser("extract reads/alignments from a single cell")
    parser.add_argument("file", help="A SAM or FASTQ file")
    parser.add_argument("barcode", help="barcode of the cell to extract")
    args = parser.parse_args()
    extract_fn = extract_barcode_sam if args.file.endswith(".sam") else extract_barcode_fastq
    extract_fn(args.file, args.barcode)
