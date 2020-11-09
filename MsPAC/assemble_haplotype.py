#!/bin/env python
import sys
import os
import pysam
import gzip
from Bio import SeqIO
from MsPAC import Pipeline

class PhasedBlocks():
    def __init__(self,line,min_phased_block,haps_to_assemble):
        self.line = line
        line = line.rstrip().split('\t')
        self.chrom = line[0]
        self.start = int(line[1])
        self.end = int(line[2])
        self.hap = str(line[3])
        self.comp_cost = line[4]
        self.length = self.end - self.start
        self.assemble = True
        if self.length < int(min_phased_block):
            self.assemble = False
        if self.hap not in haps_to_assemble.split(','):
            self.assemble = False

class HaplotypeAssembly(Pipeline):
    def __init__(self,configfile):
        Pipeline.__init__(self,configfile,"assembly")
        
    def get_read_groups_regions(self):
        read_groups = ["0","1","2","0_1","0_2"]
        regions = dict.fromkeys(read_groups,{})
        samfile = pysam.AlignmentFile(self.phased_bamfile)
        for read in samfile:
            if read.is_secondary:
                continue
            if read.is_supplementary:
                continue
            if read.is_unmapped:
                continue
            read_group = read.get_tag("RG")
            chrom = samfile.get_reference_name(read.reference_id)
            ref_start = read.reference_start
            ref_end = read.reference_end
            if read_group == "0":
                rgs = ["0","0_1","0_2"]
            if read_group == "1":
                rgs = ["1","0_1"]
            if read_group == "2":
                rgs = ["2","0_2"]
            for rg in rgs:
                if chrom not in regions[rg]:
                    regions[rg][chrom] = []
                regions[rg][chrom].append((ref_start,ref_end))
        return regions

    def merge_regions(self,regions):
        sorted_by_lower_bound = sorted(regions, key=lambda tup: tup[0])
        merged_regions = []
        for higher in sorted_by_lower_bound:
            if len(merged_regions) == 0:
                merged_regions.append(higher)
                continue
            lower = merged_regions[-1]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged_regions[-1] = (lower[0], upper_bound) 
            else:
                merged_regions.append(higher)
        return merged_regions

    def break_long_regions(self,regions):
        broken_regions = []
        for start,end in regions:
            for new_start in range(start,end,self.max_block_length):
                new_end = new_start + self.max_block_length
                if new_end > end:
                    new_end = end
                broken_regions.append((new_start,new_end))
        return broken_regions

    # Max block length
    def create_phased_bedfile(self):
        window_size = 1000000 # 10 MB //Everything should get low
        read_group_regions = self.get_read_groups_regions()
        with open(self.phased_bedfile,'w') as fh:
            for read_group in read_group_regions:
                for chrom in read_group_regions[read_group]:
                    merged_regions_pre_broken = self.merge_regions(read_group_regions[read_group][chrom])
                    if self.max_block_length == None:
                        merged_regions = merged_regions_pre_broken
                    else:
                        merged_regions = self.break_long_regions(merged_regions_pre_broken)
                    for start,end in merged_regions:
                        comp_cost = "low"
                        if end - start > window_size:
                            comp_cost = "high"
                        out = [chrom,start,end,read_group,comp_cost]
                        fh.write("%s\n" % "\t".join(map(str,out)))
                        
    def load_haplotype_blocks(self):
        if not self.non_emptyfile(self.phased_bedfile):
            self.create_phased_bedfile()
        with open(self.phased_bedfile,'r') as fh:
            for line in fh:
                phased_block = PhasedBlocks(line,self.min_phased_block,self.haps_to_assemble)
                if phased_block.assemble == False:
                    continue
                self.windows_to_assemble.append(phased_block)
                
    def assemble_windows(self):
        for window in self.windows_to_assemble:
            directory = "%s/%s/%s/%s_%s" % (self.assembly_directory,window.hap,
                                            window.chrom,window.start,window.end)
            self.create_directory(directory)
            if os.path.isfile("%s/done" % directory):
                continue
            template_bash = "%s/assemble_window_v3.sh" % self.package_bash_directory
            bashfile = "%s/assemble_window_v3.sh" % directory
            params = {
                'output': directory,
                'raw_reads_dir': self.raw_reads_directory,
                'python_scripts': self.package_python_directory,
                'hap': window.hap,
                'subreads_to_ref': self.phased_bamfile,
                'chrom': window.chrom,
                'start': window.start,
                'end': window.end,
                'threads': self.job_threads[window.comp_cost],
                'memory': int(self.job_threads[window.comp_cost])*int(self.job_memory[window.comp_cost]),
                'size': max(6000,window.end - window.start),
                'tech': self.tech
                }
            self.write_to_bashfile(template_bash,bashfile,params)
            self.jobs.append((bashfile,window.comp_cost))
        self.submitjobs()

    def merge_sequences(self):
        records = {}
        for hap in ["hap0","hap1","hap2"]:
            records[hap] = {}
            records[hap]["fa"] = []
            records[hap]["fq"] = []
        for window in self.windows_to_assemble:
            directory = "%s/%s/%s/%s_%s" % (self.assembly_directory,window.hap,
                                            window.chrom,window.start,window.end)
            raw_contigs = "%s/canu/raw.contigs.fasta" % directory
            #raw_contigs = "%s/canu/raw.quivered.contigs.fasta" % directory
            if self.non_emptyfile(raw_contigs):
                for index,record in enumerate(SeqIO.parse(raw_contigs, "fasta")):
                    record.id = "%s.%s.%s.%s.raw.%s/0/0_0" % (window.chrom,window.start,window.end,window.hap,index)
                    record.description = ""
                    record.name = ""
                    for hap,hap_key in zip(["0","1","2"],["hap0","hap1","hap2"]):
                        if hap in window.hap:
                            records[hap_key]["fa"].append(record)
            raw_contigs = "%s/canu/raw.quivered.contigs.fastq" % directory
            if self.non_emptyfile(raw_contigs):
                print raw_contigs
                for index,record in enumerate(SeqIO.parse(raw_contigs, "fastq")):
                    record.id = "%s.%s.%s.%s.raw.%s/0/0_0" % (window.chrom,window.start,window.end,window.hap,index)
                    record.description = ""
                    record.name = ""
                    for hap,hap_key in zip(["0","1","2"],["hap0","hap1","hap2"]):
                        if hap in window.hap:
                            records[hap_key]["fq"].append(record)
        SeqIO.write(records["hap0"]["fa"],self.hap0_assembly_fa,"fasta")
        SeqIO.write(records["hap1"]["fa"],self.hap1_assembly_fa,"fasta")
        SeqIO.write(records["hap2"]["fa"],self.hap2_assembly_fa,"fasta")
        SeqIO.write(records["hap0"]["fq"],self.hap0_assembly_fq,"fastq")
        SeqIO.write(records["hap1"]["fq"],self.hap1_assembly_fq,"fastq")
        SeqIO.write(records["hap2"]["fq"],self.hap2_assembly_fq,"fastq")


    def run(self):
        self.configure()
        self.load_haplotype_blocks()        
        self.assemble_windows()
        self.merge_sequences()
