def swapAttributes(gtflines, source, replace):
    # XXX not finished
    newlines = []
    repdict = dict(zip(replace, source))
    for line in gtflines:
        attrdict = attributeToDict(line)
        newdict = attrdict.copy()
        for r, s in repdict.items():
            newdict[r] = attrdict[s]
            line[r] = attrdict[s]
        line['attribute'] = buildAttributeFieldFromDict(newdict)
        newlines.append(line)
            
    return newlines

def delAttributes(gtflines, delete):
    newlines = []
    for line in gtflines:
        attrdict = attributeToDict(line)
        for x in delete:
            if x in attrdict:
                del attrdict[x]
            if x in line:
                del line[x]
            line['attribute'] = buildAttributeFieldFromDict(attrdict)
        newlines.append(line)
    return newlines

def GTFtoDict(infn):
    print "Parsing the GTF file %s." %(infn)
    gtflines = []
    infile = open(infn, 'r')
    numlines = 0
    for line in infile:
        numlines = numlines + 1
        gtflines.append(parseGTFlineToDict(line))

    print "Processed %d lines in %s." %(numlines, infn)
    infile.close()

    return gtflines

def buildAttributeFieldFromDict(attrdict):
    z = [x[0] + " " + x[1] for x in attrdict.items()]
    return("; ".join(z) + "\n")
        
def attributeToDict(linedict):
    attributes = linedict["attribute"].split(";")
    del attributes[-1]
    attributes = [x.strip() for x in attributes]

    return dict([x.strip().split(" ") for x in attributes])

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
