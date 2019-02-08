# MsPAC
**Phase reads, assemble haplotypes and detect SVs**

MsPAC takes in long reads and phased SNVs to separate the reads into two haplotypes, and assembles both haplotypes and detects structural variants. 

## Usage
MsPac is split into four steps. For each step, the input is a configuration file. A description of the configuration file is below.
#### `phase-bam`
In the first step `phase-bam`, a bam file is created. This bam file is a copy of the input bam file with a read group annotation added to the reads. A read group annotation of 1 and 2 corresponds to haplotype 1 and 2. The read group annotation of 0 corresponds to unassignable reads.
#### `prep-reads`
In the second step `prep-reads`, several bam files are created. These bam contain the raw reads seperated by chromosome and haplotype. It makes the process of searching for these reads much faster during the Quiver process, where haplotype specific reads are used to clean the haplotype-specific contigs.
#### `assembly`
In the third step `assembly`, the haplotypes are assembled. During this process folders will be created for each region. Within each folder there is a bash script that runs the assembly process. MsPAC can submit these bash scripts as a single job into the cluster (this speeds up the process).
#### `sv-calling`
In the last step `sv-calling`, the haplotypes and reference are aligned and the SVs are called. In this step, new directories will be made that holds the multiple sequence alignment and a BED file with the SVs.

```
MsPAC phase-bam run.cfg
MsPAC prep-reads run.cfg
MsPAC assembly run.cfg
MsPAC sv-calling run.cfg
```
#### Notes
1. Only works for Linux. This is due to packages only being avaible for Linux.

## Installing
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

## Running the tests
```
cd testing
sh run.sh
```
## Configuration file
All the inputs are self-explanatory except `BAM fofn` and `Phased bedfile`. `BAM fofn` is file with a list of all the raw reads in BAM format (this is the output of a PacBio run). `Phased bedfile` is a bedfile that contains 5 entries. An entry is coordinate with the haplotype and then low or high. For example: `chr1 1 100 0_1 low`. This entries are the regions that will be assembled. There is an example configuration file in the testing folder.
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
