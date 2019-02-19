# MsPAC
**Phase reads, assemble haplotypes and detect SVs**

[Introduction](#introduction)  
[Tool requirements](#tool-requirements)  
[Installation](#installation)  
[Test runs](#test-runs)<br/>
[Configuration File](#cfg-file)<br/>
[Quick Start](#quick-start)       
[Explanation of steps](#explanation-of-steps)     
[Example of output](#example-of-output)


## Introduction
MsPAC takes in long reads and phased SNVs to separate the reads into two haplotypes, and assembles both haplotypes and detects structural variants. The output is a fasta file containing both haplotypes and VCF file with SVs. The SVs are annotated with their type, size, genotype and reference, haplotype 1 and haplotype 2 sequence.

## Tool requirements
1. Linux operating system
2. [Conda package](https://conda.io/en/latest/)
3. [cluster python package](https://github.com/oscarlr/cluster)

## Installation
```
### Installing MsPAC and it's dependencies
git clone https://github.com/oscarlr/MsPAC.git
cd MsPAC
conda env create -f environment.yml 
conda activate MsPAC
python setup.py install

### Installing cluster package that's needed
cd ..
git clone https://github.com/oscarlr/cluster.git
cd cluster
python setup.py install
```

## Test runs
```
cd testing
sh run.sh
```
## Configuration File
Explanation of configuration file entries is [here](cfg_readme.md).
```
[Input]
directory = 

[Phase-bam input files]
phased vcf = 
reads aligned = 

[Phase-bam params]
sample name in VCF = 
output phased bamfile = 

[Prep reads params]
BAM fofn = 
Raw reads directory =

[Assembly params]
Minimum phased block length = 1000
Comma-seperated list of haplotypes = 0_1,0_2
Assembly directory = 
Flanking length = 1000
Phased bedfile = None

[SV calling params]
SV calling directory =
reference = 

[Other params]
cluster = No

[HIGH INTENSITY JOB]
walltime = 24
threads = 1
memory = 8
queue = private

[LOW INTENSITY JOB]
walltime = 24
threads = 1
memory = 8
queue = private
```

## Quick Start
```
MsPAC phase-bam run.cfg
MsPAC prep-reads run.cfg
MsPAC assembly run.cfg
MsPAC sv-calling run.cfg
```

## Explanation of steps
MsPac is split into four steps. For each step, the input is a configuration file. A description of the configuration file is [here](cfg_readme.md).
#### `phase-bam`
In the first step `phase-bam`, a bam file is created. This bam file is a copy of the input bam file with a read group annotation added to the reads. A read group annotation of 1 and 2 corresponds to haplotype 1 and 2. The read group annotation of 0 corresponds to unassignable reads.
#### `prep-reads`
In the second step `prep-reads`, several bam files are created. These bam contain the raw reads seperated by chromosome and haplotype. It makes the process of searching for these reads much faster during the Quiver process, where haplotype specific reads are used to clean the haplotype-specific contigs.
#### `assembly`
In the third step `assembly`, the haplotypes are assembled. During this process folders will be created for each region. Within each folder there is a bash script that runs the assembly process. MsPAC can submit these bash scripts as a single job into the cluster (this speeds up the process).
#### `sv-calling`
In the last step `sv-calling`, the haplotypes and reference are aligned and the SVs are called. In this step, new directories will be made that holds the multiple sequence alignment and a BED file with the SVs.

## Example of output
```
chr22	16610019	16610020	INS	1|0	46	46.6780821918	46.84	.	CACTGCTGTTGGGTTCTCTTTGTTTTTCCTCACAAAGGATTCCACA	.	18270	18316	/sc/orga/work/rodrio10/software/in_github/MsPAC/testing/MsPAC/sv_calling/chr22/16595201_16611082/msa.clu
```
The columns are:
```
1. chromosome
2. SV start
3. SV end
4. SV type
5. SV genotype
6. SV size
7. Haplotype 1 SV quality score 
8. Haplotype 2 SV quality score
9. Reference sequence
10. Haplotype 1 sequence
11. Haplotype 2 sequence
12. Start index position of SV in multiple sequence alignment file 
13. End index position of SV in multiple sequence alignment file 
14. Full path of multiple sequence alignment file
```
