import logging
from argparse import ArgumentParser
import os
from gtfUtils import calculateLengths, outputLengths, GTFtoDict
from gtfUtils import filterByMinLength, filterByMaxLength, outputGTF

def main():
    logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s',
                        level=logging.INFO)

    description = "Counts up size of transcripts in a GTF file. Can also " \
                  "filter on minimum transcript size."

    parser = ArgumentParser(description=description)

    parser.add_argument("-g", "--gtf", dest="gtf", default=False,
                        type=str, required=True,
                        help="gtf file to analyze")
    parser.add_argument("-l", "--length", dest="length", default=False,
                        action="store_true",
                        help="output length of each transcript")
    parser.add_argument("-m", "--min_size", dest="min_size", default=False,
                        type=int,
                        help="remove transcripts below this size.")
    parser.add_argument("-M", "--max_size", dest="max_size", default=False,
                        type=int,
                        help="remove transcripts greater than this size.")
    parser.add_argument("-o", "--outfile", dest="outfn", default=False,
                        type=str)

    args = parser.parse_args()

    if not os.path.isfile(args.gtf):
        logging.error("%s cannot be found." %(args.gtf))
        parser.print_help()
        exit(-1)

    gtflines = GTFtoDict(args.gtf)

    def checkOutFile(args):
        if not args.outfn:
            logging.error("need to provide an output filename.")
            parser.print_help()
            exit(-1)
        if os.path.isfile(args.outfn):
            logging.error("%s already exists, aborting." %(args.outfn))
            exit(-1)

    if args.length:
        checkOutFile(args)
        lengths = calculateLengths(gtflines)
        outputLengths(lengths, args.outfn)
        exit(1)

    if args.min_size:
        checkOutFile(args)
        gtflines = filterByMinLength(gtflines, args.min_size)
        
    if args.max_size:
        checkOutFile(args)
        gtflines = filterByMaxLength(gtflines, args.max_size)

    outputGTF(gtflines, args.outfn)
    
if __name__ == "__main__":
    main()
