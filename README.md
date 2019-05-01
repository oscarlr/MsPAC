# MsPAC
**Phase reads, assemble haplotypes and detect SVs**

[Introduction](#introduction)  
[Tool requirements](#tool-requirements)  
[Installation](#installation)  
[Cluster configuration](#cluster-configuration)<br/>
[Test runs](#test-runs)<br/>
[Configuration File](#configuration-file)<br/>
[Quick Start](#quick-start)       
[Explanation of steps](#explanation-of-steps)     
[Example of output](#example-of-output)<br/>
[Manuscript results](#manuscript-results)


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
## Cluster configuration
If you don't want to use the cluster use this command before running MsPAC:
```
export SJOB_DEFALLOC=""
```
If you want to use the cluster, edit the `lsf/cluster/config.py` script in `https://github.com/oscarlr/cluster.git`. The cluster package reads from this file the default configurations to run jobs in the cluster as wells as the account to use when submitting jobs. After you edit `lsf/cluster/config.py` reinstall the package using `python setup.py install` in the cluster folder.

## Test runs
```
export SJOB_DEFALLOC=""
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
### BED SV output
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
### Assembled fasta haplotype 
```
>22.16050007.16697745.0_1.raw.0/0/0_0
GACCATGTGAAACTAAGGACAACTTCAGAGCTTCACACAGCTTCAACACTGGAGAGAAAA
CAGTGAACCCACAGAAAACATCCTACAGACTGGGAGAAAATTATGGAAAACTGTGGATCT
GGAAGGGCTTCTTATCTAACATATTCAAGAAACTAATGGTCCTAAGTGGACAAAAACCAA
TATACAATGCTTGTCACACCTAAGTGGACAAAAACCAATACTAAAAATGCCCAAAAGACT
GCGTAGGCATTTCTGAAAAAACCTGAAACAGCCTCTCAGGTAACAGAAGTTTCTCCACAT
CAAGAAGAGTTTCTCCCCAGAGAACGAGTATGACCAGAAAACAGCAATAAAACTTTGGAA
TAAGAGATAAGGGCAGTGTAGATTTGCAGACAGAGGAACTATTACATACTACCTGGTTTG
AATGCAAATTTGTATACCCACTGGGAAACAGCTGGAGGTTTCTGAAACAATTAACAACAC
AACCACCAGTTCCTCTAGCCATCCCACACTGGGTATACCTGCAAAGCCAAGGAAACCTAC
```
The fasta header has the region that was assembled with the corresponding haplotype.

### Phased BAM file
```
m150131_015113_42163R_c100780292550000001823166508251570_s1_p0/8761/14473_27456
16
22
16050008
0
15543H97S27=1I1=1X12=1I20=2D13=1I8=1I9=3I6=1I11=1D5=1D1=1D11=1I3=1D11=1I4=1X1=1D3=1I8=5I10=1I1=1I3=1I5=3I15=1I8=1D38=1X16=1X5=1I8=1X4=1I12=1I4=1I14=1I1=1I8=1D6=1I4=1I13=1I15=1D1=1I21=2I12=1I1=1X6=1X3=1D5=
*
0
12343
CCAATCTCCTGGCAGCCACGCAGCCGGTCGAGAAATTTCGTCACTTGTGGCGGGTTCCCAAGCCTGTTGCCATGCAGCCTCTGGAAAGAGATCTGATTAAGTCCCAGGACTTCAGAAGAGCTGTTGCGACCTTGGCCAATGTCACTTCCTCCTTCAGGAATTGCAGTGGGCCTTAAGTGCCTTCCTCTCGGGCCCACTGGTTAT
--*%..)-.-/,.//*-/..,/.%+%"'-)(./*./)'/"(..%,(.,(+&)"*.(-,+-./-"///*/,.//+.//,/).(+)/*.//+/)//,//#/'./.,,./)/'&%///.+.////.(%/./.+.,*/+)(..$/////,+...,,/&&..((%.(/.////,$,/*'/.+//,//%-././..,'###(&(+($',,
AS:i:-54803
XS:i:15544
XE:i:28527
qs:i:15544
qe:i:28527
zm:i:-1
XL:i:12884
XT:i:1
NM:i:0
FI:i:14476
XQ:i:42999
iq:Z:0/,(121/113-223+/212121/,-#'//202*12*'2#)10%0(0.),&.#*0)0--/03.#232*3-1332223.2*111,2+123,3232.32#3'120/.23)2-3&3333.1223232%313211-*3-*(01$22332--101/.2&-01)0)0)313322-$-2*'30,23-22%/121211-'###,'/,
dq:Z:222'22*2222222222222222&2&(22+)2222222222222-2222222)2222222222222222222222222222(222222222*22222222222222222('222222222222)222222222222222222222222222222'222)'222222222222222222222222222222222(2+2)2
sq:Z:<<<<<<<<<<<<<9<<<<<<<<<<<<<<<<<<<<<<<<<%<<<<<<<<<<<<<<<<<<<<<<<$<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<5<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<3<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
mq:Z:03/6>0>75,834;=1280-*:3,C.-1/*:@9C3@1}>}}19}:857}:}&}:753:564:C}>9M323=62'4>1;:}94''A?=9:}5C;6}?=}=}FG6-71;B4}DI5>6-)=D<3@17}9;40'7}D:1-}95};57=8};/A9-9B}}8:346@(8283C@7}<8}}88}963=6}74;5A8>:7&%-.}9<
st:Z:AACCGAGAAGTTACTAACATACTAATTGATCTCCCGGGATGACAGGTGTTATTTGGAAACCTAAGTGGTAACGTACTAAGAGTTCCCTCTCGAGTCGGCCTGAAACTTCAGGACTCCTCTAGTGGTATCAAGGTTAACCGTGACAGGAAGAAGGACTTCCGGTACTGTTTAAGGCCTGTAAGGAAGAGATTTAAACAGT
dt:Z:NNNANNCNNNNNNNNNNNNNNNNANTTNNTANNNNNNNNNNNNNTNNNNNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNNANNNNNNNNNANNNNNNNNNNNNNNNNNGTNNNNNNNNNNNNGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNNNCCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNANGN
ip:Z:S26,94,16,11,14,14,77,19,28,28,26,91,24,5,35,8,11,46,10,31,75,27,16,257,26,54,44,15,39,47,16,7,37,7,17,57,50,20,26,30,23,49,63,6,636,54,12,33,33,16,15,141,360,18,14,25,42,8,17,21,51,10,17,34,19,13,12
RG:Z:2
```
RG tag has the haplotype information. This read can be accessed using the `samtools view` command with the `-r` option. For example, `samtools view reads.bam 22:16050008-16050108 -r 2`.
## Manuscript results
[HG002 haplotype assemblies](https://rodrio10.u.hpc.mssm.edu/MsPAC/hg002_assembly/haplotypes/)
