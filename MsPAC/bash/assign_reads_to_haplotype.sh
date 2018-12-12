#!/bin/bash
set -e -x

python_scripts=$1
vcffile=$2
bamfile=$3
phased_bamfile=$4
vcf_sample_name=$5
##

if [ ! -s "${phased_bamfile}.bai" ]
then
    # python ${python_scripts}/assign_reads_to_haplotypes.py \
    # 	${vcffile} \
    # 	${bamfile} \
    # 	${vcf_sample_name} \
    # 	${phased_bamfile}
    samtools index ${phased_bamfile}
fi


