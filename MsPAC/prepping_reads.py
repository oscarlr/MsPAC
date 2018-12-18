#!/bin/env python
import sys
import pysam
from MsPAC import Pipeline


class PrepReads(Pipeline):
    def __init__(self,configfile):
        Pipeline.__init__(self,configfile,"prep-reads")

    def get_reads_per_group(self):
        print "\tGetting reads per chrom per group..."
        reads_per_read_group = {}
        samfile = pysam.AlignmentFile(self.phased_bamfile)
        reads_mapped = samfile.mapped
        count = 0
        for read in samfile:            
            count += 1.0
            read_group = read.get_tag("RG")
            chrom = samfile.get_reference_name(read.reference_id)
            if chrom not in reads_per_read_group:
                reads_per_read_group[chrom] = {}
            if read_group not in reads_per_read_group[chrom]:
                reads_per_read_group[chrom][read_group] = []
            reads_per_read_group[chrom][read_group].append(read.query_name)            
            if count % 50000  == 0:
                print "\t\tStatus: %s" % (count/reads_mapped)
        return reads_per_read_group

    def extract_raw_reads(self,read_names,inbamfile):
        raw_reads = []
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
                    raw_reads.append(x)
        return raw_reads

    def get_raw_reads(self,read_names):        
        headers = {}
        raw_reads = []
        with open(self.raw_reads_in_bam_format,'r') as fofh:
            for i,fn in enumerate(fofh):
                print "\t\tLooking into bam file #: %s" % i
                fn = fn.rstrip()
                inbamfile = pysam.AlignmentFile(fn,check_sq=False)
                reads = self.extract_raw_reads(read_names,inbamfile)
                if len(reads) > 0:
                    headers[fn] = dict(inbamfile.header)
                inbamfile.close()
                for read in reads:
                    raw_reads.append(read)
        return raw_reads,headers

    def split_raw_reads_into_groups(self,raw_reads,max_num_reads_per_file,header,outdir):
        for index,i in enumerate(range(0,len(raw_reads),max_num_reads_per_file)):
            reads = raw_reads[i:i + max_num_reads_per_file]
            outbamfilefn = "%s/%s.bam" % (outdir,index)
            outbamfile = pysam.AlignmentFile(outbamfilefn,"wb", header=header)
            for read in reads:
                outbamfile.write(read)
            outbamfile.close()
            pysam.index(outbamfilefn)

    def merge_headers(self,headers):
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

    def split_raw_reads(self):
        max_num_reads = 200000 # 200,000 reads is around 5GB
        read_names_per_read_group = self.get_reads_per_group()
        for chrom in read_names_per_read_group:
            for read_group in read_names_per_read_group[chrom]:        
                print "\tWorking on chrom %s and read group %s" % (chrom,read_group)
                print "\tExtracting %s reads" % len(read_names_per_read_group[chrom][read_group])
                outdir = "%s/%s/%s" % (self.raw_reads_directory,chrom,read_group)
                self.create_directory(outdir)
                read_names = read_names_per_read_group[chrom][read_group]
                raw_reads,headers = self.get_raw_reads(read_names)
                header = self.merge_headers(headers)
                self.split_raw_reads_into_groups(raw_reads,max_num_reads,header,outdir)

    def run(self):
        self.configure()
        self.split_raw_reads()
