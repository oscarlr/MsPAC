[Example](#example)  
[Explanation](#explanation)

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
Technology = ONT

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
`sample name in VCF` is the sample name in the VCF file. `output phased bamfile` is the output BAM file with the input reads phased. In the `output phased bamfile`, the input reads have a new read group tag. The read group tag `0`, `1`, and `2` correspond to unphased reads, haplotype 1, and haplotype 2.

```
[Prep reads params]
BAM fofn = reads.fofn
Raw reads directory = MsPAC/prep_reads
```
`BAM fofn` is a file created by MsPAC that list all the BAM files created by the `prep-reads` step. `Raw reads directory` is the directory with the BAM files `prep-reads` creates. The BAM files contain the raw reads seperated by chromosome.

```
[Assembly params]
Minimum phased block length = 1000
Comma-seperated list of haplotypes = 0_1,0_2
Assembly directory = MsPAC/assembly
Flanking length = 1000
Phased bedfile = None
```
`Minimum phased block length` is the minimum size that will be assembled. `Comma-seperated list of haplotypes` are the haplotypes that will be assembled. The options are: `0`,`1`,`2`,`0_1`, and `0_2`. `0`,`1`, and `2` are unambiguous regions, haplotype 1 and haplotype 2. `0_1` and `0_2` are haplotype 1 and 2 with the reads from unambiguous regions added to both haplotype 1 and 2. `Assembly directory ` is the directory with the regions assembled. `Flanking length` is an extra amount of bases added to both ends of each region to be assembled. `Phased bedfile` is a bed file with the regions to assemble. It is created by MsPAC if none is given. `Phased bedfile` should have this tab-delimited format:
`chromosome start end haplotype low/high`, for example:
```
22	16050007	16697745	1	low
22	16847850	17262375	1	low
22	17262464	18711525	1	low
22	18712024	18712281	1	low
22	50414777	51244565	0	low
22	16050007	16697745	2	low
22	50414777	51244565	2	low
22	16050007	16697745	0_2	low
22	20609570	50364777	0_1	high
22	50414777	51244565	0_1	low
```

```
[SV calling params]
SV calling directory = MsPAC/sv_calling
reference = input_data/chr22.fa
```
`SV calling directory` is the directory with the output from the `sv-calling` MsPAC step. `reference` is the genome reference in fasta format. 

```
[Other params]
cluster = No
```
If `cluster` is "Yes", then the assembly and sv-calling jobs will be sent to the cluster.

```
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
Regions labelled `low` in `Phased bedfile` will use the `[LOW INTENSITY JOB]` configuration, and similarly for `high` regions.
