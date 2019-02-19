[Example](#example)  
[Explanation]

# Example
```
[Input]
directory = MsPAC

[Phase-bam input files]
phased vcf = input_data/test.vcf.gz
reads aligned = input_data/test.bam

[Phase-bam params]
sample name in VCF = 20977
output phased bamfile = input_data/test_phased.bam

[Prep reads params]
BAM fofn = reads.fofn
Raw reads directory = MsPAC/prep_reads

[Assembly params]
Minimum phased block length = 1000
Comma-seperated list of haplotypes = 0_1,0_2
Assembly directory = MsPAC/assembly
Flanking length = 1000
Phased bedfile = None

[SV calling params]
SV calling directory = MsPAC/sv_calling
reference = input_data/chr22.fa

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

# Explanation
```
[Input]
directory = MsPAC
```
The directory that MsPAC writes to.

```
[Phase-bam input files]
phased vcf = input_data/test.vcf.gz
reads aligned = input_data/test.bam
```
`phase vcf` is the input phased VCF file with the phased SNPs. `reads aligned` is the input BAM file with reads aligned. 

```
[Phase-bam params]
sample name in VCF = 20977
output phased bamfile = input_data/test_phased.bam
```
`sample name in VCF` is the sample name in the VCF file. `output phased bamfile` is the output BAM file with the input reads phased. The input reads have a new read group tag. The read group tag `0` is unphased reads. The read group tag `1` is haplotype 1 and the tag `2` is haplotype 2.
