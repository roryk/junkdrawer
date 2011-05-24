#!/bin/bash
if [ $# -lt 2 ]; then
    SCRIPT_NAME=`basename $0`
    echo -e "usage: $SCRIPT_NAME gtffile outfile \n"
    echo "gtffile must have the gene_id attribute and the transcript_id
next to each other, like this: gene_id foo; transcript_id bar; "
    exit
fi

echo "Counting up the number of isoforms for each gene..."
egrep -o -e 'gene_id "\w+"; transcript_id "\w+"' $1 | cut -f2,4 -d" " | tr -d ";\"" | uniq | cut -f1 -d" " | uniq -c | awk '{print $1}' | sort | uniq -c | sort -k2,2n > $2

