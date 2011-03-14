# extract sequences for splice junctions mapped with tophat from the
# .juncs file
from optparse import OptionParser
from subprocess import call
import time
import os

def parseTophatJunctions(infn, outfn):
    infile = open(infn, 'r')
    outfile = open(outfn, 'w')
    junctions = 0
    for line in infile:
        if 'track name' not in line:
            junctions = junctions + 1
            bl = line.split('\t')
            blocksizes = bl[10].split(',')
            blockstarts = bl[11].split(',')
            blockstarts = [int(x) for x in blockstarts]
            blocksizes = [int(x) for x in blocksizes]
            bl[1] = int(bl[1])
            block1start = bl[1] + blockstarts[0]
            block1end = bl[1] + blockstarts[0] + blocksizes[0]
            block2start = bl[1] + blockstarts[1]
            block2end = bl[1] + blockstarts[1] + blocksizes[1]
            outline = (bl[0] + "\t" + str(block1start) + "\t" + str(block1end) +
                       "\t" + bl[3] + "\t" + bl[4] + "\t" + bl[5])
            outfile.write(outline + '\n')
            outline = (bl[0] + "\t" + str(block2start) + "\t" + str(block2end) +
                       "\t" + bl[3] + "\t" + bl[4] + "\t" + bl[5])
            outfile.write(outline + '\n')
    infile.close()
    outfile.close()
    print "Read in %d junctions." %(junctions)

def combineJunctionEnds(tmpfilename, outprefix):
    """
    combines reads with the same junction name to put together the
    whole splice junction sequence
    """
    tmpfilein = open(tmpfilename, 'r')
    juncdict = {}
    seqlen = 0
    junctions = 0

    for line in tmpfilein:
        if ">" in line:
            jname = line.split(">", 2)[1]
        else:
            seqlen = seqlen + len(line.rstrip())
            if jname not in juncdict:
                junctions = junctions + 1
                juncdict[jname] = line.rstrip()
            else:
                juncdict[jname] = juncdict[jname] + line.rstrip()

    tmpfilein.close()
    outfile = open(outprefix + ".fa", 'w')
    for key in juncdict:
        header = ">" + key + "\n"
        seq = juncdict[key] + "\n"
        outfile.write(header)
        outfile.write(seq)
    outfile.close()

    print "Wrote %d junctions of %d total bases." %(junctions, seqlen)

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", dest="filename", 
                      help="tophat junctions bed file")
    parser.add_option("-o", dest="outfile", 
                       help="output suffix")
    parser.add_option("-g", dest="genomefile", 
                       help="genome filename")
    (options, args) = parser.parse_args()

    if options.filename is None:
        print parser.print_help()
        exit(1)
    if not os.path.isfile(options.filename):
        print "%s cannot be found." %(options.filename)
    if options.outfile is None:
        print parser.print_help()
        exit(1)
    if options.genomefile is None:
        print parser.print_help()
        exit(1)
    if not os.path.isfile(options.genomefile):
        print "%s cannot be found." %(options.genomefile)
        
        
    tmpfilejuncs = "tmpjunctionfile" + str(time.time())
    parseTophatJunctions(options.filename, tmpfilejuncs)
    tmpfileseq = "tmpseqfile" + str(time.time())
                                           
    print "Looking up sequences in %s with fastaFromBed." %(options.genomefile)
    call("fastaFromBed -fi %s -bed %s -s -name -fo %s" %(options.genomefile, tmpfilejuncs, tmpfileseq), shell=True)

    print "Combining junction ends."
    combineJunctionEnds(tmpfileseq, options.outfile)

    print "Cleaning up temporary files."
    os.remove(tmpfilejuncs)

if __name__ == "__main__":
    main()
