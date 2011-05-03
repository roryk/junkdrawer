from gtfUtils import GTFtoDict, outputGTF, swapAttributes, delAttributes
from gtfUtils import addAttribute, filterAttributes
from argparse import ArgumentParser
import os
import logging

def main():

    logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s',
                        level=logging.INFO)

    description = "Replaces the attributes in -r with the ones from -s. " \
                  "Deletes the attributes in -d."
    parser = ArgumentParser(description=description)
    parser.add_argument("-g", "--gtf", dest="gtf", default=False,
                        type=str, required=True,
                        help="combined.gtf file from Cuffcompare")
    parser.add_argument("-s", "--source", nargs="*", dest="source",
                        default=False,
                        help="attributes to replace the -r attributes")
    parser.add_argument("-r", "--replace", nargs="*", dest="replace",
                        default=False,
                        help="attributes to be replaced")
    parser.add_argument("-d", "--delete", nargs="*", dest="delete",
                        default=False,
                        help="attributes to be deleted")
    parser.add_argument("-a", "--add", dest="add", default=False,
                        type=str,
                        help="file of attributes to be added")
    parser.add_argument("-f", "--filter", dest="filter", default=False,
                        type=str,
                        help="file of attribute to filter on")
    parser.add_argument("-o", "--output", dest="output", default=False,
                        type=str, required=True,
                        help="output filename")
    args = parser.parse_args()

    if not os.path.isfile(args.gtf):
        logging.error("%s cannot be found." %(args.gtf))
        parser.print_help()
        exit(-1)

    if os.path.isfile(args.output):
        logging.error("%s already exists." %(args.output))
        parser.print_help()
        exit(-1)

    # check to make sure arguments make sense
    if args.source or args.replace:
        if len(args.source) != len(args.replace):
            logging.error("Source and replacement lengths must be the same.")
            parser.print_help()
            exit(-1)

    if not (args.source or args.replace or args.delete or
            args.add or args.filter):
        logging.error("Must provide at least one action.")
        parser.print_help()
        exit(-1)

    gtflines = GTFtoDict(args.gtf)
    if args.source:
        logging.info("Swapping attributes.")
        gtflines = swapAttributes(gtflines, args.source, args.replace)

    if args.add:
        if not os.path.isfile(args.add):
            logging.error("%s cannot be found." %(args.add))
            exit(-1)
        logging.info("Adding attributes.")
        gtflines = addAttribute(gtflines, args.add)
        
    if args.delete:
        logging.info("Deleting attributes.")
        gtflines = delAttributes(gtflines, args.delete)

    if args.filter:
        logging.info("Filtering by attribute in file %s." %(args.filter))
        gtflines = filterAttributes(gtflines, args.filter)
        
    outputGTF(gtflines, args.output)
    
if __name__ == "__main__":
    main()
