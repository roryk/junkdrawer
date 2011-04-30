from gtfUtils import GTFtoDict, outputGTF, swapAttributes, deleteAttributes
from argparse import ArgumentParser

def main():
    description = "In a GTF file with attributes, replaces the attributes "\
                  "specified by the -r argument with the arguments from " \
                  "-s. These are specified by a space separated list. " \
                  " Attributes in the -d option are deleted from " \
                  " the GTF file."

    parser = ArgumentParser(description)
    parser.add_argument("-f", nargs=1, required=True, metavar="gtf",
                        help="combined.gtf file from Cuffcompare")
    parser.add_argument("-s", nargs="*", metavar="source", default=False,
                        help="attributes to replace the -r attributes")
    parser.add_argument("-r", nargs="*", metavar="replace", default=False,
                        help="attributes to be replaced")
    parser.add_argument("-d", nargs="*", metavar="delete", default=False,
                        help="attributes to be deleted")

    args = parser.parse_args()

    # check to make sure arguments make sense
    if args.source or args.replace:
        if len(args.source) != len(args.replace):
            print "Source and replacement lengths must be the same."
            parser.print_help()
            exit(-1)

    if not args.source or args.replace or args.delete:
        print "Must provide at least one action."
        parser.print_help()
        exit(-1)

    gtflines = GTFtoDict(args.gtf)
    if args.source:
        print "Swapping attributes."
        gtflines = swapAttributes(gtflines, args.source, args.replace)
        
    if args.delete:
        print "Deleting atributes."
        gtflines = deleteAttributes(gtflines, args.delete)

    outputGTF(gtflines)
    
if __name__ == "__main__":
    main()
