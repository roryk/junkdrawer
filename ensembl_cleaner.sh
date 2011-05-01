#!/bin/bash
# This script takes a GTF file from Ensembl and cleans it up for use with
# Tophat and Cufflinks.

if [ $# -lt 2 ]; then
    echo "usage: $0 ensembl.gtf outfile_cleaned.gtf"
    exit 1
fi

echo "Prefacing chromosome names with chr..."
cat $1 | awk -F"\t" '{OFS="\t"; $1="chr"$1;print}' > $1.tmp
echo "Replacing chrMT with chrM..."
cat $1.tmp | awk  -F"\t" '{OFS="\t"; if($1=="chrMT") $1="chrM"; print}' > $2

echo "Cleaning up temporary files..."
rm $1.tmp 
