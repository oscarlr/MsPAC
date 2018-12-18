#!/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import math
import pysam
import os

from MsPAC import Pipeline

class SVCaller(Pipeline):
    def __init__(self,configfile):
        Pipeline.__init__(self,configfile,"sv-calling")

    def create_new_record(self,record,i,start,end):     
        new_record_name = "%s.%s" % (i,record.id)
        new_record_seq = record.seq[start:end]
        new_record_qual = record.letter_annotations["phred_quality"][start:end]
        new_record = SeqRecord(new_record_seq,id=new_record_name,name=new_record_name,description="")
        new_record.letter_annotations["phred_quality"] = new_record_qual
        return new_record
    
    def create_new_fastq(self,infastq,outfastq,window):
        records = []            
        for record in SeqIO.parse(infastq, "fastq"):
            if len(record.seq) > window:
                num_windows = int(math.ceil(len(record.seq)/window))
                size_of_window = len(record.seq)/num_windows
                i = 0
                for i in range(num_windows - 1):
                    start = i*size_of_window
                    end = (i+1)*size_of_window
                    new_record = create_new_record(record,i,start,end)
                    records.append(new_record)
                start = (i+1)*size_of_window
                new_record = create_new_record(record,i+1,start,None)
                records.append(new_record)
            else:
                records.append(record)
        SeqIO.write(records,outfastq,"fastq")        

    def split_fastq_files(self):
        window = 500000
        for hap in ["hap1","hap2"]:
            if hap == "hap1":
                infastq = self.hap1_assembly_fq
                outfastq = self.hap1_assembly_split_fq
            elif hap == "hap2":
                infastq = self.hap2_assembly_fq
                outfastq = self.hap2_assembly_split_fq
            self.create_new_fastq(infastq,outfastq,window)

    def map_assembly(self):
        template_bash = "%s/map.sh" % self.package_bash_directory
        for hap in ["hap1","hap2"]:
            if hap == "hap1":
                fastq = self.hap1_assembly_split_fq
            elif hap == "hap2":
                fastq = self.hap2_assembly_split_fq
            prefix = "%s/%s_to_ref" % (self.sv_calling_directory,hap)
            bashfile = "%s/map_%s_assembly.sh" % (self.sv_calling_directory,hap)
            params = {
                'input': fastq,
                'ref': self.reference,
                'prefix': prefix,
                'threads': self.job_threads["high"]
                }
            self.write_to_bashfile(template_bash,bashfile,params)
            self.jobs.append((bashfile,self.job_threads["high"]))
        self.run_locally()

    def get_msa_coordinates(self):
        if os.path.isfile("%s/msa_coords.bed" % self.sv_calling_directory):
            return
        template_bash = "%s/get_msa_coords.sh" % self.package_bash_directory
        bashfile = "%s/get_msa_coords.sh" % self.sv_calling_directory
        params = {
            'python_packages': self.package_python_directory,
            'ref': self.reference,
            'sv_calling_dir': self.sv_calling_directory
            }
        self.write_to_bashfile(template_bash,bashfile,params)
        self.jobs.append((bashfile,self.job_threads["high"]))
        self.run_locally()
    
    def read_bedfile(self,bedfile):
        regions = {}
        with open(bedfile,'r') as bedfh:
            for line in bedfh:
                line = line.rstrip().split('\t')
                chrom = line[0]
                coord = (int(line[1]),int(line[2]))
                if chrom not in regions:
                    regions[chrom] = []
                regions[chrom].append(coord)
        return regions

    def get_hap_sequence(self,hap_bamfile,regions):
        samfile = pysam.AlignmentFile(hap_bamfile,'rb')
        hap_sequences = {}
        for chrom in regions:
            for start,end in regions[chrom]:
                for contig in samfile.fetch(chrom,start,end):
                    if contig.is_unmapped:
                        continue
                    if contig.is_secondary:
                        continue
                    if contig.is_supplementary:
                        continue
                    if contig.reference_start > start:
                        continue
                    if contig.reference_end < end:
                        continue
                    aligned_pairs = contig.get_aligned_pairs()
                    query_start = None
                    query_end = None
                    matched_ref_start = None
                    matched_ref_end = None
                    for query_pos, ref_pos in aligned_pairs:
                        if query_pos == None:
                            continue
                        if ref_pos == None:
                            continue
                        if int(ref_pos) <= int(start):
                            query_start = query_pos
                            matched_ref_start = ref_pos
                        query_end = query_pos
                        matched_ref_end = ref_pos
                        if int(ref_pos) > int(end):
                            break
                    assert query_start != None
                    assert query_end != None
                    hap_sequence = contig.query_sequence[query_start:query_end]
                    if contig.query_qualities != None:
                        sequence_qual = contig.query_qualities[query_start:query_end]
                    else:
                        sequence_qual = [0]*len(hap_sequence)
                    assert len(sequence_qual) == len(hap_sequence)
                    hap_sequences[(chrom,start,end)] = (hap_sequence,sequence_qual)
            return hap_sequences

    def extract_msa_sequence(self):
        msa_coordsfn = "%s/msa_coords.bed" % self.sv_calling_directory
        msa_coords = self.read_bedfile(msa_coordsfn)
        hap1_bamfn = "%s/hap1_to_ref.sorted.bam" % self.sv_calling_directory
        hap1_sequence = self.get_hap_sequence(hap1_bamfn,msa_coords)
        hap2_bamfn = "%s/hap2_to_ref.sorted.bam" % self.sv_calling_directory
        hap2_sequence = self.get_hap_sequence(hap2_bamfn,msa_coords)
        fasta = pysam.FastaFile(self.reference)
        for chrom in msa_coords:
            for i,(start,end) in enumerate(msa_coords[chrom]):
                ref_seq = fasta.fetch(reference=chrom,start=max(1,start),end=end) 
                h1_seq,h1_qual = hap1_sequence[(chrom,start,end)]
                h2_seq,h2_qual = hap2_sequence[(chrom,start,end)]
                directory = "%s/%s/%s_%s" % (self.sv_calling_directory,chrom,start,end)
                self.create_directory(directory)
                outseqfn = "%s/seq.fa" % directory
                outqualfn = "%s/seq.qual" % directory
                with open(outseqfn,'w') as outfh:
                    outfh.write(">ref\n%s\n" % ref_seq)
                    outfh.write(">hap1\n%s\n" % h1_seq)
                    outfh.write(">hap2\n%s\n" % h2_seq)
                with open(outqualfn,'w') as outqualfh:
                    outqualfh.write(">hap1\n%s\n" % ",".join(map(str,h1_qual)))
                    outqualfh.write(">hap2\n%s\n" % ",".join(map(str,h2_qual)))

    def calls_svs_from_msa(self):
        msa_coordsfn = "%s/msa_coords.bed" % self.sv_calling_directory
        msa_coords = self.read_bedfile(msa_coordsfn)
        template_bash = "%s/sv_calling.sh" % self.package_bash_directory
        for chrom in msa_coords:
            for i,(start,end) in enumerate(msa_coords[chrom]):
                directory = "%s/%s/%s_%s" % (self.sv_calling_directory,chrom,start,end)
                bashfile = "%s/sv_calling.sh" % directory
                params = {
                    'dir': directory,
                    'python_scripts': self.package_python_directory,
                    'chrom': chrom,
                    'start': start,
                    'end': end
                    }
                self.write_to_bashfile(template_bash,bashfile,params)
                self.jobs.append((bashfile,self.job_threads["low"]))
        self.submitjobs()        

    def run(self):
        self.configure()
        self.split_fastq_files()
        self.map_assembly()
        self.get_msa_coordinates()
        self.extract_msa_sequence()
        self.calls_svs_from_msa()
