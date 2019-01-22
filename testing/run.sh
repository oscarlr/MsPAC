#!/bin/bash
set -e -x

if [ -z "${SJOB_DEFALLOC}" ]
then
    export SJOB_DEFALLOC=""
fi

mkdir -p input_data
cd input_data

if [ ! -s reads.bam.bai ]
then
    curl -O https://rodrio10.u.hpc.mssm.edu/MsPAC/reads.bam
    curl -O https://rodrio10.u.hpc.mssm.edu/MsPAC/reads.bam.bai
fi

if [ ! -s test.bam.bai ]
then
    curl -O https://rodrio10.u.hpc.mssm.edu/MsPAC/test.bam
    curl -O https://rodrio10.u.hpc.mssm.edu/MsPAC/test.bam.bai
fi

if [ ! -s test.vcf.gz.tbi ]
then
    curl -O https://rodrio10.u.hpc.mssm.edu/MsPAC/test.vcf.gz
    curl -O https://rodrio10.u.hpc.mssm.edu/MsPAC/test.vcf.gz.tbi
fi

if [ ! -s chr22.fa ]
then
    curl -O http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr22.fa.gz
    gunzip chr22.fa.gz
    samtools faidx chr22.fa
fi
cd -

ls input_data/reads.bam > reads.fofn

MsPAC phase-bam run.cfg
MsPAC prep-reads run.cfg
MsPAC assembly run.cfg
MsPAC sv-calling run.cfg
