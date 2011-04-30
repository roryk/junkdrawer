def swapAttributes(gtflines, source, replace):
    # XXX not finished
    return gtflines

def delAttributes(gtflines, delete):
    # XXX not finished
    return gtflines

def GTFtoDict(infn):
    print "Parsing the GTF file %s."
    gtflines = []
    infile = open(infn, 'r')
    numlines = 0
    for line in infile:
        numlines = numlines + 1
        gtflines.append(parseGTFlineToDict(line))

    print "Processed %d lines in %s." %(numlines, infn)
    infile.close()

    return gtflines

def addAttributesToGTFline(linedict):
    attributes = linedict["attribute"].split(";")
    del attributes[-1]
    attributes = [x.strip() for x in attributes]
    
    for attrib in attributes:
        attr = attrib.strip().split(" ")
        linedict[attr[0]] = attr[1].strip("\"")
    return linedict

def parseGTFlineToDict(line):
    values = line.split("\t")
    keys = ["seqname", "source", "feature", "start", "end", "score",
            "strand", "unknown", "attribute"]
    linedict = dict(zip(keys, values))
    linedict = addAttributesToGTFline(linedict)
    
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
                             line['attribute']])
        outfile.write(outline)
        written = written + 1

    outfile.close()
    print "Wrote %d lines to %s." %(written, outfn)
