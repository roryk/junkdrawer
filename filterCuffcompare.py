from optparse import OptionParser
from sets import Set
import os

def outputTracking(tracklines, outfn):
    print "Writing tracking file."
    outfile = open(outfn, 'w')
    written = 0
    samples = [str(x) + "_samp" for x in range(0, tracklines[0]['nsamples'])]
    for line in tracklines:
        samplist = [line[x] for x in samples]
        outline = "\t".join([line['transid'], line['locusid'],
                             line['ref_geneid'], line['ref_transid']] +
                            samplist) + "\n"
        outfile.write(outline)
        written = written + 1
    outfile.close()
    print "Wrote %d lines to %s." %(written, outfn)

def parseTracklineToDict(line):
    line = line.split("\t")
    # first 4 columns are mandatory, the rest contain the sample info
    nsamples = len(line) - 4
    line = [nsamples] + line
    samples = [str(x) + "_samp" for x in range(0, nsamples)]
    keys = ["nsamples", "transid", "locusid", "ref_geneid",
            "ref_transid"] + samples
    
    linedict = dict(zip(keys, line))
    return linedict

def outputGTF(gtflines, outfn):
    print "Writing GTF file."
    outfile = open(outfn, 'w')
    written = 0
    for line in gtflines:
        
        outline = "\t".join([line['seqname'], line['source'],
                             line['feature'], line['start'],
                             line['end'], line['score'],
                             line['strand'], line['unknown'],
                             line['attribute']]) + "\n"
        outfile.write(outline)
        written = written + 1

    outfile.close()
    print "Wrote %d lines to %s." %(written, outfn)
    
def parseGTFlineToDict(line):
    values = line.split("\t")
    keys = ["seqname", "source", "feature", "start", "end", "score",
            "strand", "unknown", "attribute"]
    linedict = dict(zip(keys, values))
    return linedict

def parseTophatGTFlineToDict(line):
    linedict = parseGTFlineToDict(line)
    attributes = linedict["attribute"].split(";")
    del attributes[-1]
    attributes = [x.strip() for x in attributes]
    
    for attrib in attributes:
        attr = attrib.strip().split(" ")
        linedict[attr[0]] = attr[1].strip("\"")
    return linedict

def parseTophatGTF(gtf):
    print "Parsing GTF file."
    gtflines = []
    gfile = open(gtf, 'r')
    numlines = 0
    for line in gfile:
        numlines = numlines + 1
        gtflines.append(parseTophatGTFlineToDict(line))

    print "Processed %d lines in %s." %(numlines, gtf)
    gfile.close()

    return gtflines
    
def parseTracking(track):
    print "Parsing tracking file."
    tracklines = []
    tfile = open(track, 'r')
    numlines = 0
    for line in tfile:
        numlines = numlines + 1
        tracklines.append(parseTracklineToDict(line))

    print "Processed %d lines in %s." %(numlines, track)
    tfile.close()

    return tracklines
         
def filterGTFByTracking(gtflines, tracklines):
    idset = Set()
    kept = 0
    processed = 0
    for trackline in tracklines:
        idset.add(trackline['transid'])        

    newgtf = []
    for line in gtflines:
        processed = processed + 1
        if line['transcript_id'] in idset:
            newgtf.append(line)
            kept = kept + 1
            
    print "Kept %d out of %d lines in the GTF file." %(kept, processed)
    return(newgtf)

def filterByNumsamples(tracklines, threshold):
    
    print "Filtering for lines with at least %d samples of " \
          "support." %(threshold)
    keptlines = []
    kept = 0
    for linedict in tracklines:
        if enoughSamples(linedict, threshold):
            kept = kept + 1
            keptlines.append(linedict)

    print "%d lines had at least %d samples of support." %(kept, threshold)
    return keptlines


def enoughSamples(linedict, threshold):
    """
    returns true if the number of samples supporting the transcript is
    greater than or equal to the threshold and false otherwise
    """
    support = linedict.values().count("-")
    return(support >= threshold)
        
               
def main():
    parser = OptionParser()
    parser.add_option("-f", dest="gtf", default=None,
                      help="combined.gtf file from Cuffcompare")
    parser.add_option("-t", dest="tracking", default=None,
                      help="tracking file from Cuffcompare")
    parser.add_option("-s", dest="samples", default=None,
                      help="number of samples a ID must appear in")

    (options, args) = parser.parse_args()
    
    if options.gtf is None:
        print parser.print_help()
        exit(-1)
    if options.tracking is None:
        print parser.print_help()
        exit(-1)
    if options.samples is None:
        print parser.print_help()
        exit(-1)

    trackout = options.tracking + ".filtered"
    gtfout = options.gtf + ".filtered"

    if os.path.isfile(trackout):
        print "%s already exists, aborting." %(trackout)
        exit(-1)
    if os.path.isfile(gtfout):
        print "%s already exists, aborting." %(trackout)
        exit(-1)
    
    tracklines = parseTracking(options.tracking)
    tracklines = filterByNumsamples(tracklines, int(options.samples))
    gtflines = parseTophatGTF(options.gtf)
    gtflines = filterGTFByTracking(gtflines, tracklines)
    outputTracking(tracklines, trackout)
    outputGTF(gtflines, gtfout)
if __name__ == "__main__":
    main()
    
