# extract sequences for splice junctions mapped with tophat from the
# .juncs file
from optparse import OptionParser
import tophatUtils as th
from subprocess import call
import os

def parseExtendedJuncLineToDict(line):
    linedict = {}
    bl = line.split('\t')
    idline = bl[4].split(':')
    linedict['chrom'] = idline[0]
    linedict['left'] = idline[1]
    linedict['right'] = idline[2]
    linedict['strand'] = idline[3]
    linedict['score'] = bl[5]
    return linedict

def formatExtendedJunc(linedict, offset):
    """
    returns two lines, one for each side of the junction
    """
    chrom = linedict['chrom']
    left = str(linedict['left'])
    right = str(linedict['right'])
    left_off = str(linedict['left'] - offset)
    right_off = str(linedict['right'] + offset)
    strand = linedict['strand']
    score = str(linedict['score'])

    id = ":".join([chrom, left, right, strand])
    outleft = "\t".join([chrom, left_off, left, strand, id, score])
    outright = "\t".join([chrom, right, right_off, strand, id, score])
    return (outleft, outright)

def convertTophatBedToExtendedJuncs(infn, outfn, offset):
    infile = open(infn, 'r')
    outfile = open(outfn, 'w')
    junctions = 0

    for line in infile:
        if 'track name' not in line:
            junctions = junctions + 1
            linedict = th.parseTophatBedLineToDict(line)
            linedict = th.calculateJunctionEnds(linedict)
            (outleft, outright) = formatExtendedJunc(linedict, offset)
            outfile.write(outleft + '\n')
            outfile.write(outright + '\n')
            
    infile.close()
    outfile.close()
    print "Processed %d junctions." %(junctions)

def convertExtendedJuncsToSeq(infn, outfn, genfn):

    call("fastaFromBed -fi %s -bed %s -name -fo %s" %(infn, outfn, genfn))

    
def combineJunctionEnds(seqfn):
    """
    combines reads with the same junction name to put together the
    whole splice junction sequence
    """
    seqfile = open(seqfn, 'r')
    juncdict = {}
    seqlen = 0
    junctions = 0

    for line in seqfile:
        if ">" in line:
            jname = line.split(">", 2)[1]
        else:
            seqlen = seqlen + len(line.rstrip())
            if jname not in juncdict:
                junctions = junctions + 1
                juncdict[jname] = line.rstrip()
            else:
                juncdict[jname] = juncdict[jname] + line.rstrip()

    seqfile.close()
    seqfile = open(seqfn, 'w')
    for key in juncdict:
        header = ">" + key 
        seq = juncdict[key] 
        seqfile.write(header)
        seqfile.write(seq)
        seqfile.write("\n")
    seqfile.close()

    print "Wrote %d junctions of %d total bases." %(junctions, seqlen)

def main():
    #usage = "usage: python %prog [options]"
    parser = OptionParser()
    parser.add_option("-f", dest="infn", default=None,
                      help="tophat junctions bed file")
    parser.add_option("-g", dest="genomefile", 
                       help="genome filename (FASTA format)")
    parser.add_option("-t", dest="keeplen", default=0,
                       help="keep only this many bases per side of junction")

    (options, args) = parser.parse_args()
    
    if options.infn is None:
        print parser.print_help()
        exit(-1)
    if not os.path.isfile(options.infn):
        print "%s cannot be found." %(options.infn)
        exit(-1)
    if options.genomefile is None:
        if not os.path.isfile(options.genomefile):
            print "%s cannot be found." %(options.genomefile)
            exit(-1)

    if options.keeplen is 0:
        print parser.print_help()

    juncfn = options.infn + ".juncs"
    seqfn = options.infn + ".fa"

    if os.path.isfile(juncfn):
        print "%s already exists, aborting." %(juncfn)

    if os.path.isfile(seqfn):
        print "%s already exists, aborting." %(seqfn)

    print "Reading in junctions."
    convertTophatBedToExtendedJuncs(options.infn, juncfn, options.keeplen)
    print "Looking up sequences in %s with fastaFromBed." %(options.genomefile)
    convertExtendedJuncsToSeq(juncfn, seqfn, options.genomefile)
    print "Combining the junction ends."
    combineJunctionEnds(seqfn)
    
if __name__ == "__main__":
    main()
