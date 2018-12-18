#!/bin/bash
set -e -x

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

cd -

MsPAC phase-bam run.cfg
MsPAC prep-reads run.cfg