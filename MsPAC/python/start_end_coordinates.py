#!/bin/env python
import pysam
import sys

bamfile = sys.argv[1]
samfile = pysam.AlignmentFile(bamfile)
for i,read in enumerate(samfile):
    if read.is_unmapped:
        continue
    if read.is_secondary:
        continue
    if read.is_supplementary:
        continue
    output = [samfile.getrname(read.reference_id),read.reference_start,read.reference_end - 1,read.query_name]
    print "\t".join(map(str,output))
        
