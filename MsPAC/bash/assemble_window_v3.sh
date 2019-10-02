#!/bin/bash
set -x

output=$1
bash=$2
smrtsuite=$3
threads=$4
bams_fofn=$5
CANU=$6
size=$7
subreads_to_ref=$8
hap=$9
reference=${10}
memory=${11}
raw_reads_dir=${12}
start=${13}
end=${14}
chrom=${15}
python_scripts=${16}

if [ "${hap}" == "0_1" ]
then
    samtools view ${subreads_to_ref} -r 0 ${chrom}:${start}-${end} | awk '{ print ">"$1"\n"$10}' > ${output}/reads.fasta
    samtools view ${subreads_to_ref} -r 1 ${chrom}:${start}-${end} | awk '{ print ">"$1"\n"$10}' >> ${output}/reads.fasta
elif [ "${hap}" == "0_2" ]
then
    samtools view ${subreads_to_ref} -r 0 ${chrom}:${start}-${end} | awk '{ print ">"$1"\n"$10}' > ${output}/reads.fasta
    samtools view ${subreads_to_ref} -r 2 ${chrom}:${start}-${end} | awk '{ print ">"$1"\n"$10}' >> ${output}/reads.fasta
else
    samtools view ${subreads_to_ref} -r ${hap} ${chrom}:${start}-${end} | awk '{ print ">"$1"\n"$10}' > ${output}/reads.fasta
fi

samtools faidx ${output}/reads.fasta

if [ ! -s ${output}/canu/raw.contigs.fasta ]
then    
    ##rm -fr ${output}/canu
    canu \
	-p raw \
	-d ${output}/canu \
	contigFilter="2 1000 1.0 1.0 2" \
	corMinCoverage=0 \
	stopOnLowCoverage=0 \
	minThreads=${threads} \
	genomeSize=${size} \
	useGrid=0 \
	-pacbio-raw ${output}/reads.fasta
fi

if [  -s ${output}/canu/raw.contigs.fasta ]
then
    samtools faidx ${output}/canu/raw.contigs.fasta
    if [ ! -s ${output}/subreads.bam ]
    then
	if [ "${hap}" == "0_1" ]
	then
	    ls ${raw_reads_dir}/${chrom}/1/*bam > ${output}/subreads.fofn
	    ls ${raw_reads_dir}/${chrom}/0/*bam >> ${output}/subreads.fofn
	elif [ "${hap}" == "0_2" ]
	then
	    ls ${raw_reads_dir}/${chrom}/2/*bam > ${output}/subreads.fofn
	    ls ${raw_reads_dir}/${chrom}/0/*bam >> ${output}/subreads.fofn
	else
	    ls ${raw_reads_dir}/${chrom}/${hap}/*bam > ${output}/subreads.fofn
	fi
	if [ ! -s ${output}/subreads.fofn ]
	then
	    echo "" > ${output}/done
	    exit 0
	fi
	cut -f1 ${output}/reads.fasta.fai > ${output}/reads.id
	python \
	    ${python_scripts}/extract_raw_reads_from_bam_fofn.py \
	    ${output}/reads.id \
	    ${output}/subreads.fofn \
	    ${output}/subreads.bam
	pbindex ${output}/subreads.bam
    fi
    if [ ! -s ${output}/canu/reads_to_canu_contigs.sorted.bam.pbi ]
    then
        blasr \
	    ${output}/subreads.bam \
            ${output}/canu/raw.contigs.fasta \
            --bestn 1 \
            --bam \
            --nproc ${threads} \
            --out ${output}/canu/reads_to_canu_contigs.bam
        samtools sort -@ ${threads} ${output}/canu/reads_to_canu_contigs.bam -o ${output}/canu/reads_to_canu_contigs.sorted.bam
        pbindex ${output}/canu/reads_to_canu_contigs.sorted.bam
    fi
    if [ ! -s ${output}/canu/raw.quivered.contigs.fastq ]
    then
        samtools faidx ${output}/canu/raw.contigs.fasta
        arrow \
            --referenceFilename ${output}/canu/raw.contigs.fasta \
            -j ${threads} \
            -o ${output}/canu/raw.quivered.contigs.fastq \
            -o ${output}/canu/raw.quivered.contigs.fasta \
            ${output}/canu/reads_to_canu_contigs.sorted.bam
    fi
fi

echo "" > ${output}/done