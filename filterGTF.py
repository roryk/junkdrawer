import sys
from sets import Set
import logging
import os
from argparse import ArgumentParser
from gtfUtils import GTFtoDict, formatGTFLine, addAttributesToGTFline

def filterAttributes(gtflines, B, F, v=False):
    Bfile = open(B, 'r')
    filter_set = Set()
    for line in Bfile:
        filter_set.add(line.strip())

    matched = []
    notmatched = []

    for line in gtflines:
        linedict = addAttributesToGTFline(line)
        if linedict[F] in filter_set:
            matched.append(line)
        else:
            notmatched.append(line)

    if v:
        return notmatched
    else:
        return matched
        

def main():

    logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s.',
                        level=logging.INFO)

    description = "Keeps GTF lines of file A that have values in field F in file B." \
                  "The -v option removes the lines instead of keeping them."
    parser = ArgumentParser(description=description)
    parser.add_argument('A', metavar='A',
                        help='GTF file to be filtered')
    parser.add_argument('F', metavar='F',
                        help='field to filter on')
    parser.add_argument('B', metavar='B',
                        help='file with filter items')
    parser.add_argument('-v', dest='inverse', default=False,
                        action='store_true',
                        help='exclude lines instead of keeping')

    args = parser.parse_args()

    if not os.path.isfile(args.A):
        logging.error("%s cannot be found." %(args.A))
        parser.print_help()
        exit(-1)

    if not os.path.isfile(args.B):
        logging.error("%s cannot be found." %(args.B))
        parser.print_help()
        exit(-1)

    gtflines = GTFtoDict(args.A)
    gtflines = filterAttributes(gtflines, args.B, args.F, args.inverse)

    logging.info("Writing GTF file.")
    written = 0
    for line in gtflines:
        outline = formatGTFLine(line)
        sys.stdout.write(outline)
        written = written + 1
        
    logging.info("Wrote %d lines." %(written))
    
if __name__ == "__main__":
    main()
