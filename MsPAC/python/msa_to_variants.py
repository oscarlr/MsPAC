#!/bin/env python
#from MsPAC.python_scripts.hmm import *
from hmm2 import *
from pomegranate import *
from Bio import AlignIO
import numpy as np
import pysam
import sys

def get_msa_sequence(clufn):
    alignment = AlignIO.read(clufn,"clustal")
    sequences = {}
    for i,sequence in enumerate(alignment):        
        if sequence.id not in ["ref","hap1","hap2"]:
            if i == 0:
                sequence.id = "ref"
            if i == 1:
                sequence.id = "hap1"
            if i == 2:
                sequence.id = "hap2"
        sequences[sequence.id] = str(sequence.seq).upper()
    return sequences

def get_observations(sequence):
    obs = []
    if len(sequence) == 3:
        h1_seq = sequence["hap1"]
        h2_seq = sequence["hap2"]
        chrom_seq = sequence["ref"]
        for h1,h2,chrom in zip(h1_seq,h2_seq,chrom_seq):
            obs.append(observations["3"][h1][h2][chrom])
        return obs

def get_quality_scores(sequence,quality_scores,start_index,end_index,padding):
    quality_start = start_index - sequence[:start_index].count("-") #- padding
    quality_end = end_index - sequence[:end_index].count("-") #+ padding
    #print quality_start,quality_end
    if quality_end - quality_start < (padding*2):
        quality_start = max(0,quality_start - padding)
        quality_end = quality_end + padding
    scores = quality_scores[quality_start:quality_end]
    #print scores,start_index,end_index,quality_start,quality_end
    mean = float(sum(scores))/len(scores)
    return mean

def three_sequence_msa_variants(sequence,path,clufn,chrom,ref_start,ref_end,quality_scores,padding):
    start = False
    current_sv = None
    current_sv_start_index = None
    current_sv_current_index = None
    current_sv_hap1_seq = []    
    current_sv_hap2_seq = []
    current_sv_ref_seq = []
    ref_index = -1
    ref_sv_start = None
    ref_sv_index = None
    for i,(h1,h2,ref,p) in enumerate(zip(sequence["hap1"],sequence["hap2"],sequence["ref"],path[1:-1])):
        state_index, state = p
        if ref != "-":
            ref_index += 1
        if start == False:
            if "-" not in (h1,h2,ref,p) and state.name == "NORMAL":
                start = True
            continue
        if state.name != "NORMAL":
            if current_sv == None:
                current_sv = state.name
                current_sv_start_index = i
                current_sv_current_index = i
                ref_sv_start = ref_index
                ref_sv_index = ref_index
                continue
            current_sv_current_index += 1 
            ref_sv_index = ref_index
            if h1 != "-":
                current_sv_hap1_seq.append(h1)
            if h2 != "-":
                current_sv_hap2_seq.append(h2)
            if ref != "-":
                current_sv_ref_seq.append(ref)
        if i != current_sv_current_index and current_sv != None:
            if i + 1 == len(sequence["hap1"]) and "-" in (h1,h2,ref):
                continue
            sv_type, genotype = current_sv.split("_")
            sv_len = max(len(current_sv_hap1_seq),len(current_sv_hap2_seq),len(current_sv_ref_seq))
            if quality_scores != None:
                hap1_qual_score = get_quality_scores(sequence["hap1"],quality_scores["hap1"],current_sv_start_index,current_sv_current_index,padding)
                hap2_qual_score = get_quality_scores(sequence["hap2"],quality_scores["hap2"],current_sv_start_index,current_sv_current_index,padding)
            else:
                hap1_qual_score = 60
                hap2_qual_score = 60
            output = [chrom,                              
                      ref_sv_start + ref_start, # 0-based/UCSC Genome format
                      ref_sv_index + ref_start + 1, # 0-based/UCSC Genome format    
                      sv_type,
                      genotype,
                      sv_len,
                      hap1_qual_score,
                      hap2_qual_score,
                      "".join(current_sv_ref_seq) if len("".join(current_sv_ref_seq)) > 0 else ".",
                      "".join(current_sv_hap1_seq) if len("".join(current_sv_hap1_seq)) > 0 else ".",
                      "".join(current_sv_hap2_seq) if len("".join(current_sv_hap2_seq)) > 0 else ".",
                      current_sv_start_index,
                      current_sv_current_index,
                      clufn] 
            print "\t".join(map(str,output))
            current_sv = None
            current_sv_hap1_seq = []    
            current_sv_hap2_seq = []
            current_sv_ref_seq = []

def path_to_variants(path,sequence,clufn,chrom,start,end,quality_scores,padding):
    if len(sequence) == 3:
        three_sequence_msa_variants(sequence,path,clufn,chrom,start,end,quality_scores,padding)

def load_qual_scores(qual_scoresfn):
    qual_scores = {}
    with open(qual_scoresfn,'r') as fh:
        for line in fh:
            if ">" in line:
                name = line[1:].rstrip()
            else:
                qual =  line.rstrip().split(',')
                qual_scores[name] = map(int,qual)
    return qual_scores

def main():
    clufn = sys.argv[1]
    chrom = sys.argv[2]
    start = int(sys.argv[3])
    end = int(sys.argv[4])
    if sys.argv[5] != "None":
        qual_scoresfn = sys.argv[5]
        qual_scores = load_qual_scores(qual_scoresfn)
    else:
        qual_scores = None
    if sys.argv[6] != "None":
        padding = int(sys.argv[6])
    else:
        padding = 5
    sequence = get_msa_sequence(clufn)
    obs = get_observations(sequence)
    log, path = model[str(len(sequence))].viterbi(obs)
    path_to_variants(path,sequence,clufn,chrom,start,end,qual_scores,padding)   

if __name__ == "__main__":
    main()
