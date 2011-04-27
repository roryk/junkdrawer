# converts a Tophat junctions.bed file to juncs format
from optparse import OptionParser
import tophatUtils as th

def formatJuncs(linedict):
    """
    outputs a line in the juncs format
    """
    chrom = linedict['chrom']
    left = str(linedict['left'])
    right = str(linedict['right'])
    strand = linedict['strand']
    score = str(linedict['score'])

    id = ":".join([chrom, left, right, strand])
    outline = "\t".join([chrom, left, right, strand, id, score])

    return outline

def convertTophatBedToJuncs(infn, outfn):
    infile = open(infn, 'r')
    outfile = open(outfn, 'w')
    junctions = 0
    for line in infile:
        if 'track name' not in line:
            junctions = junctions + 1
            linedict = th.parseTophatBedLineToDict(line)
            linedict = th.calculateJunctionEnds(linedict)
            outline = formatJuncs(linedict)
            outfile.write(outline)

    infile.close()
    outfile.close()
    print "Processed %d junctions." %(junctions)

def main():
    #usage = "usage: python %prog [options]"
    parser = OptionParser()
    parser.add_option("-f", dest="infn", default=None,
                      help="tophat junctions bed file")
    parser.add_option("-o", dest="outfn", 
                       help="output junctions file")

    (options, args) = parser.parse_args()
    
    if options.infn is None:
        print parser.print_help()
        exit(-1)
    if options.outfn is None:
        print parser.print_help()
        exit(-1)

    convertTophatBedToJuncs(options.infn, options.outfn)

if __name__ == "__main__":
    main()
