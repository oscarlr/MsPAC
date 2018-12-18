#!/bin/bash
set -x

if [ ! -s ${prefix}.sorted.bam.bai ]
then
    blasr \
	${input} \
	${ref} \
	--bestn 1 \
	--bam \
	--nproc ${threads} \
	--out ${prefix}.bam
    
    samtools \
	sort -@ ${threads} \
	${prefix}.bam \
	-o ${prefix}.sorted.bam
    
    samtools index ${prefix}.sorted.bam
    
    rm -f ${prefix}.bam
fi
