# extract sequences for splice junctions mapped with tophat from the
# .juncs file
from optparse import OptionParser
from subprocess import call
import time
import os

def parseTophatBedLineToDict(line):
    linedict = {}
    bl = line.split('\t')
    blocksizes = bl[10].split(',')
    blockstarts = bl[11].split(',')
    linedict['blockstarts'] = [int(x) for x in blockstarts]
    linedict['blocksizes'] = [int(x) for x in blocksizes]
    linedict['chromStart'] = int(bl[1])
    linedict['chromEnd'] = int(bl[2])
    linedict['name'] = bl[3]
    linedict['score'] = bl[4]
    linedict['strand'] = bl[5]
    linedict['chrom'] = bl[0]
    return(linedict)

def convertTophatBedToJuncs(infn, outfn):
    infile = open(infn, 'r')
    outfile = open(outfn, 'w')
    junctions = 0
    for line in infile:
        if 'track name' not in line:
            junctions = junctions + 1
            linedict = parseTophatBedLineToDict(line)
            chrom = linedict['chrom']
            left = linedict['chromStart'] + linedict['blocksizes'][0]
            right = linedict['chromEnd'] - linedict['blocksizes'][1]
            strand = linedict['strand']
            outline = (chrom + "\t" + str(left) + "\t" + str(right) + "\t"
                       + strand)
            outfile.write(outline + '\n')

    infile.close()
    outfile.close()
    print "Processed in %d junctions." %(junctions)    
            
    
def parseTophatJunctions(infn, outfn, keeplen):
    infile = open(infn, 'r')
    outfile = open(outfn, 'w')
    junctions = 0
    kept = 0
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
            # trim the reads if the option to do that has been set
            if keeplen:
                block1start = block1end - keeplen
                block2end = block2start + keeplen
            if (block1start > block1end or block2start > block2end):
                continue
            kept = kept + 1
            outline = (bl[0] + "\t" + str(block1start) + "\t" + str(block1end) +
                       "\t" + bl[3] + "\t" + bl[4] + "\t" + bl[5])
            outfile.write(outline + '\n')
            outline = (bl[0] + "\t" + str(block2start) + "\t" + str(block2end) +
                       "\t" + bl[3] + "\t" + bl[4] + "\t" + bl[5])
            outfile.write(outline + '\n')
    infile.close()
    outfile.close()
    print "Read in %d junctions." %(junctions)
    print "Kept %d junctions. (%f)" %(kept, float(kept)/float(junctions))

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
        header = ">" + key 
        seq = juncdict[key] 
        outfile.write(header)
        outfile.write(seq)
        outfile.write("\n")
    outfile.close()

    print "Wrote %d junctions of %d total bases." %(junctions, seqlen)

def main():
    #usage = "usage: python %prog [options]"
    parser = OptionParser()
    parser.add_option("-f", dest="filename", default=None,
                      help="tophat junctions bed file")
    parser.add_option("-o", dest="outfile", 
                       help="output prefix")
    parser.add_option("-g", dest="genomefile", 
                       help="genome filename (FASTA format)")
    parser.add_option("-t", dest="keeplen", default=0,
                       help="keep only this many bases per side of junction")
    (options, args) = parser.parse_args()

    if options.filename is None:
        print parser.print_help()
        exit(-1)
    if not os.path.isfile(options.filename):
        print "%s cannot be found." %(options.filename)
        exit(-1)
    if options.outfile is None:
        print parser.print_help()
        exit(-1)
    if options.genomefile is None:
        print parser.print_help()
        exit(-1)
    if not os.path.isfile(options.genomefile):
        print "%s cannot be found." %(options.genomefile)
        exit(-1)
        
        
    tmpfilejuncs = "tmpjunctionfile" + str(time.time())
    parseTophatJunctions(options.filename, tmpfilejuncs, int(options.keeplen))
    tmpfileseq = "tmpseqfile" + str(time.time())
                                           
    print "Looking up sequences in %s with fastaFromBed." %(options.genomefile)
    call("fastaFromBed -fi %s -bed %s -name -fo %s" %(options.genomefile, tmpfilejuncs, tmpfileseq), shell=True)

    print "Combining junction ends."
    combineJunctionEnds(tmpfileseq, options.outfile)

    print "Cleaning up temporary files."
    os.remove(tmpfilejuncs)
    os.remove(tmpfileseq)

if __name__ == "__main__":
    main()
