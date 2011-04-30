from gtfUtils import GTFtoDict, outputGTF, swapAttributes, delAttributes
from argparse import ArgumentParser
import os

def main():

    description = "Replaces the attributes in -r with the ones from -s. " \
                  "Deletes the attributes in -d."
    parser = ArgumentParser(description=description)
    parser.add_argument("-f", dest="gtf", default=False, type=str,
                        help="combined.gtf file from Cuffcompare")
    parser.add_argument("-s", nargs="*", dest="source", default=False,
                        help="attributes to replace the -r attributes")
    parser.add_argument("-r", nargs="*", dest="replace", default=False,
                        help="attributes to be replaced")
    parser.add_argument("-d", nargs="*", dest="delete", default=False,
                        help="attributes to be deleted")

    args = parser.parse_args()

    if args.gtf:
        if not os.path.isfile(args.gtf):
            print "%s cannot be found." %(args.gtf)
            parser.print_help()
            exit(-1)
            
    # check to make sure arguments make sense
    if args.source or args.replace:
        if len(args.source) != len(args.replace):
            print "Source and replacement lengths must be the same."
            parser.print_help()
            exit(-1)

    if not (args.source or args.replace or args.delete):
        print "Must provide at least one action."
        parser.print_help()
        exit(-1)

    gtflines = GTFtoDict(args.gtf)
    if args.source:
        print "Swapping attributes."
        gtflines = swapAttributes(gtflines, args.source, args.replace)
        
    if args.delete:
        print "Deleting attributes."
        gtflines = delAttributes(gtflines, args.delete)

    outputGTF(gtflines, args.gtf + ".edited")
    
if __name__ == "__main__":
    main()
