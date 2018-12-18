#!/bin/bash
set -e -x

if [ ! -s ${ref}.fai ]
then
    samtools faidx ${ref}    
fi

cut -f1,2 ${ref}.fai > ${sv_calling_dir}/chrom.sizes

for i in 1 2
do
    python ${python_packages}/start_end_coordinates.py \
	${sv_calling_dir}/hap${i}_to_ref.sorted.bam > ${sv_calling_dir}/hap${i}_to_ref.bed
    
    bedtools genomecov \
	-bg \
	-i ${sv_calling_dir}/hap${i}_to_ref.bed \
	-g ${sv_calling_dir}/chrom.sizes \
	| awk '$4 == 1' \
	> ${sv_calling_dir}/hap${i}_to_ref.no_overlap.bed
    
    bedtools intersect \
	-a ${sv_calling_dir}/hap${i}_to_ref.bed \
	-b ${sv_calling_dir}/hap${i}_to_ref.no_overlap.bed \
	> ${sv_calling_dir}/hap${i}_to_ref.no_overlap.contig.bed
done

bedtools intersect \
    -a ${sv_calling_dir}/hap1_to_ref.no_overlap.contig.bed \
    -b ${sv_calling_dir}/hap2_to_ref.no_overlap.contig.bed \
    | cut -f-3 \
    > ${sv_calling_dir}/msa_coords.bed \

rm -f ${sv_calling_dir}/hap1_to_ref.no_overlap.bed
rm -f ${sv_calling_dir}/hap2_to_ref.no_overlap.bed

