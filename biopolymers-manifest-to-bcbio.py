"""
script to convert manifest files generated from the Biopolymers facility into
CSV files suitable for running with bcbio-nextgen

we'll need to read the CSV file in and for each sample generate what
we expect the filename to be. then we need to dump out
filename,samplename columns for each sample, where samples that
need to be joined have the same samplename
we just need to spit out the R1 filename because bcbio-nextgen will
find the other one by walking the directories

we might have to fix bcbio-nextgen to be able to use the full path of files

/net/hsphfs1/srv/export/hsphfs1/share_root/chb/projects/felipe-rnaseq/april-round/FC_02156/Unaligned_1234_PF_mm1/Data/Project_fglehn

this has everything in it, so we can just use this as the main directory and not
worry about handling the directory names

files look like this:

LIB020500_GEN00052261_S9_L001_R2.fastq.bz2

the number 9 comes from the manifest, so we need to capture that too.

FC_ID,Run Folder, Lane, Library Pool ID, Library Pool Name, Library Id, Library Name, Index, Index Name
FC_02156,160428_NB500917_0148_AH3G53AFXX,1,LIB020500,Pool 2,GEN00052253,HC-2.1,TAAGGCGA-GCGTAAGA,N701-S517

so we can find the files by doing
{LIB_POOL_ID}_{LIBRARY_ID}_S{LINENUMBER}_L{LANE-LEFT-PADDED}_R{READ}.fastq.bz

and then each line is {FILE_NAME},{LIBRARY_NAME}

and we should be good to go
"""
from argparse import ArgumentParser


COL_LOOKUP = {"LIB_POOL_ID": 3, "LIBRARY_ID": 5, "LANE": 2, "SAMPLENAME": 6}
FILE_FORMAT = ("{LIB_POOL_ID}_{LIBRARY_ID}_S{LINENUMBER}_L{LANE_PADDED}_"
               "R1.fastq.bz")

def file_from_manifest_line(line, samplenum):
    tokens = line.split(",")
    LIB_POOL_ID = tokens[COL_LOOKUP["LIB_POOL_ID"]]
    LIBRARY_ID = tokens[COL_LOOKUP["LIBRARY_ID"]]
    LANE_PADDED = tokens[COL_LOOKUP["LANE"]].zfill(3)
    SAMPLENAME = tokens[COL_LOOKUP["SAMPLENAME"]]
    LINENUMBER = samplenum
    return FILE_FORMAT.format(**locals())

def get_samplename(line):
    tokens = line.split(",")
    return tokens[COL_LOOKUP["SAMPLENAME"]]

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("manifest", help="CSV manifest from biopolymers")
    args = parser.parse_args()
    sample_lookup = {}
    print "filename,description"
    with open(args.manifest) as in_handle:
        # skip header
        in_handle.next()
        for (linenum, line) in enumerate(in_handle):
            samplename = get_samplename(line)
            if samplename not in sample_lookup:
                sample_lookup[samplename] = linenum + 1
            samplenum = sample_lookup[samplename]
            filename = file_from_manifest_line(line, samplenum)
            print ",".join([filename, samplename])
