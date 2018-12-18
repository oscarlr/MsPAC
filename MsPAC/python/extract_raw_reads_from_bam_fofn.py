#!/bin/env python
import sys
import pysam

readsfn = sys.argv[1]
inbamfofn = sys.argv[2] 
outbamfile = sys.argv[3]

def merge_headers(headers):
    header = {}
    for fn in headers:
        if len(header) == 0:
            header["HD"] = headers[fn]["HD"]
            header["RG"] = []
            header["PG"] = []
            header["SQ"] = []
        for rg in headers[fn]["RG"]:
            header["RG"].append(rg)
        for pg in headers[fn]["PG"]:
            header["PG"].append(pg)
    return header

read_names = set()
with open(readsfn,'r') as fh:
    for line in fh:
        name = line.rstrip()        
        read_names.add(name)

outreads = []
headers = {}
with open(inbamfofn,'r') as fh:
    for i,fn in enumerate(fh):
        print ">>>>>>>>> %s" % i
        fn = fn.rstrip()
        inbamfile = pysam.AlignmentFile(fn,check_sq=False)
        read_names_indexed = pysam.IndexedReads(inbamfile)
        read_names_indexed.build()
        for name in read_names:
            try:
                read_names_indexed.find(name)
            except KeyError:
                pass
            else:
                iterator = read_names_indexed.find(name)
                for x in iterator:
                    outreads.append(x)
                    headers[fn] = dict(inbamfile.header)
        inbamfile.close()

header = merge_headers(headers)
outbamfile = pysam.AlignmentFile(outbamfile,"wb", header=header)
for read in outreads:
    outbamfile.write(read)
outbamfile.close()

