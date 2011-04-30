# converts a Tophat junctions.bed file to juncs format
from optparse import OptionParser
import tophatUtils as th
import os

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
    outline = outline + "\n"
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

    (options, args) = parser.parse_args()
    
    if options.infn is None:
        print parser.print_help()
        exit(-1)
    outfn = options.infn.rsplit(".", 1)[0] + ".juncs"
    if os.path.isfile(outfn):
        print "%s already exists." %(outfn)
        exit(-1)

    print "Converting to juncs."
    convertTophatBedToJuncs(options.infn, outfn)

if __name__ == "__main__":
    main()
