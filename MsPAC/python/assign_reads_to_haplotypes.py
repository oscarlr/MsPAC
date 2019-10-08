#!/bin/env python
import sys
import vcf
import pysam
import numpy as np

unphased_tag = 0
haplotype1_tag = 1
haplotype2_tag = 2

def read_phased_snps(vcffile,sample):
    phased_snps = {}
    vcf_reader = vcf.Reader(open(vcffile, 'r'))
    for record in vcf_reader:
        if record.var_subtype == "deletion":
            continue
        if record.var_subtype == "insertion":
            continue
        phased_snps.setdefault(record.CHROM, {})
        try:
            alleles = record.genotype(sample)['PGT'].split("|")
        except:
            alleles = record.genotype(sample)['GT'].split("|")
        if len(alleles) == 2 and alleles[0] != alleles[1]:            
            alleles = map(int, alleles)
            if max(alleles) != 1:
                continue
            allele_bases = [record.REF,record.ALT[0]]
            sample_bases = map(lambda x: allele_bases[x], alleles)
            phased_snps[record.CHROM][record.POS-1] = sample_bases
    return phased_snps

def create_tag(hap):
    haptag = ("RG", str(hap), "Z")
    return haptag

def calculate_prob(basetups,hap_dict,threshold=0.9,lr_threshold=10):
    '''
    Returns a 0, 1 or -1 as the value of a read
    '''
    log_lr_threshold = np.log(lr_threshold)
    prob0  = 1.
    prob1  = 1.
    lprob0 = 0.
    lprob1 = 0.
    for basetup in basetups:
        rpos, b, qual = basetup
        #try:
        if len(hap_dict[rpos]) == 0:
            print rpos
        b0, b1 = hap_dict[rpos]
        #except:
        #    print rpos
        #    sys.exit()
        if b == b1 or b == b0: # check to make sure the base is called
            if b == b0:
                prob0 *= (1.-10.**(-qual/10.))
                prob1 *= (10.**(-qual/10.))
                lprob0 += np.log(1.-10.**(-qual/10.))
                lprob1 += np.log(10.**(-qual/10.))
            else:
                prob1 *= (1.-10.**(-qual/10.))
                prob0 *= (10.**(-qual/10.))
                lprob1 += np.log(1.-10.**(-qual/10.))
                lprob0 += np.log(10.**(-qual/10.))
    prob01 = prob0 + prob1
    if (prob0 > 0 and prob1 > 0) and prob0 / (prob01) >= threshold:
        if lprob0 - lprob1 >= log_lr_threshold:
            return create_tag(haplotype1_tag)
    if (prob0 > 0 and prob1 > 0) and prob1 / (prob01) >= threshold:
        if lprob1 - lprob0 >= log_lr_threshold:
            return create_tag(haplotype2_tag)
    if (prob0 == 0 or prob1 == 0) and lprob0 - lprob1 >= log_lr_threshold:
        return create_tag(haplotype1_tag)
    if (prob0 == 0 or prob1 == 0) and lprob1 - lprob0 >= log_lr_threshold:
        return create_tag(haplotype2_tag)
    return create_tag(unphased_tag)

def phase_read(read,phased_snps,chrom):
    unphased_tag = 0
    haplotype1_tag = 1
    haplotype2_tag = 2
    base_tuples = []
    if chrom in phased_snps:
        for read_pos, ref_pos in read.get_aligned_pairs():
            if ref_pos in phased_snps[chrom]:
                if read_pos is not None:
                    if read.query_qualities == None:
                        qual = 8
                    else:
                        qual = read.query_qualities[read_pos]
                    base = read.query_sequence[read_pos]
                    base_tuples.append((ref_pos,base,qual))
    if len(base_tuples) > 0:
        read_group_tag = calculate_prob(base_tuples,phased_snps[chrom])
    else:
        read_group_tag = create_tag(unphased_tag)
    #print "%s\t%s\t%s" % (len(base_tuples),read.query_name,read_group_tag[1])
    read_tags = read.get_tags()
    tags_to_add = []
    for tag in read_tags:
        if tag[0] != "RG":
            tags_to_add.append(tag)
    tags_to_add.append(read_group_tag)
    read.set_tags(tags_to_add)
    return read

def create_header(unphased_bam):    
    readgroup_unphased = { "ID": unphased_tag } 
    readgroup_hap1 = { "ID": haplotype1_tag }
    readgroup_hap2 = { "ID": haplotype2_tag }
    if "RG" in unphased_bam.header:
        for group in unphased_bam.header["RG"]:
            for item in group:
                if item != "ID":
                    readgroup_unphased[item] = group[item]
                    readgroup_hap1[item] = group[item]
                    readgroup_hap2[item] = group[item]
    phased_bam_header = unphased_bam.header.to_dict()
    phased_bam_header["RG"] = [readgroup_unphased,readgroup_hap1,readgroup_hap2]
    return phased_bam_header

def main(vcffile,bamfile,vcf_sample_name,outbamfile):
    #outbamfile = "%s.phased.bam" % bamfile[0:-4]
    phased_snps = read_phased_snps(vcffile,vcf_sample_name)
    unphased_bam = pysam.AlignmentFile(bamfile,'rb')
    phased_bam_header = create_header(unphased_bam)
    phased_bam = pysam.AlignmentFile(outbamfile,'wb',header=phased_bam_header)
    for read in unphased_bam.fetch():
        if read.is_secondary:
            continue
        if read.is_unmapped:
            continue
        if read.is_supplementary:
            continue
        if not read.is_secondary:
            tagged_read = phase_read(read,phased_snps,unphased_bam.get_reference_name(read.reference_id))        
            phased_bam.write(tagged_read)
    unphased_bam.close()
    phased_bam.close()

if __name__ == '__main__':
    vcffile = sys.argv[1]
    bamfile = sys.argv[2]
    vcf_sample_name = sys.argv[3]
    phased_bamfile = sys.argv[4]
    sys.exit(main(vcffile,bamfile,vcf_sample_name,phased_bamfile))
