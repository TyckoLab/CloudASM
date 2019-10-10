
# CloudASM: a cloud-based, ultra-efficient pipeline for mapping allele-specific DNA methylation

Last updated: October 9, 2019. Please check our preprint on [biorxiv](https://www.biorxiv.org/).

## Table of contents

Here, the TOC

***********

## Overview

CloudASM is a turnkey pipeline designed to call allele-specific CpG methylation in whole methylomes. It is designed to run on [Google Cloud Platform](https://cloud.google.com/) (GCP). 

Because it leverages the Google's severless data warehouse, all steps specific to ASM-calling (as opposed to alignment, variant calling, and methylation calling), it is 500x more efficient in terms of CPU-hours when compared to a traditional server-based approach using a cluster.

This pipeline starts from the zipped fastq files hosted on a bucket and outputs a table on BigQuery of all regions in a sample showing allele-specific methylation.

## Biology significance of CloudASM

Important stuff
- we filter out the reads where the confidence in the SNP letter is less than 30
- we remove CpGs from the context file where the C or G overlap with the SNP
- we remove all SNPs that are not within 500 bp of at least a CpG


DMR: 
3 significant ASM CpG in the same direction
2 consecutive significant CpGs in the same direction
at least 20% difference between the REF reads and the ALT reads and FDR < 0.05

Bis-SNP reports SNPs in positive strand of the reference genome (it's NOT bisulfite-converted)
Bismark reports in positive strand but it is bisulfite-converted, requiring careful hand



From a biology standpoint, ASM can be calculated at the single CpG level (ASM CpG) or across a region (ASM DMR). Our pipeline uses criteria on both to output what we consider biologicaly-relevant ASM. All biological parameters relevant to ASM can be customized through the following variables:
- `GENOME`: the reference genome to be used in alignment. It can be either `hg19` or `GRCh38`.
- `DMR_EFFECT`: the minimum effect size in terms of methylation percentage between allele A and allele B of a variant.
- `CPG_COV`: the minimum coverage of each CpG in a given allele. 
- `CPG_PER_DMR`: the minimum number of CpGs per allele and the minimum number of CpGs per variant with ASM in the same direction.
- `CONSECUTIVE_CPG`: the minimum number of CpGs with significant ASM in the same direction, in a given ASM DMR.
- `SNP_SCORE`: the minimum score in ASCII required for SNP nucleotides in a read. If the score is too low, we could not tell which allele is included in the read.

The pipeline is designed to run the following steps sequentially:

1. Create a bisulfite-converted genome from the reference genome that was chosen.
1. Unzip fastq files and split them into smaller files.
1. Trim each pair of fastq files
1. Align each pair of fastq files
1. Merge BAM files per chromosome and prepare them for variant calling.
1. Net methylation calling
1. 1 

## Pipeline overview

The pipeline follows these steps:

1. 
2. Unzip fastq files and trim them in 12M-row fastq files.
3. Trim and align each pair of fastq files. Split the output BAM file in chromosome-specific BAM files.
4. Merge all BAM files per chromosome. Remove duplicates. Perform net methylation.
5. Perform SNP calling.
6. Compute allele-specific methylation.

## Bioinformatics packages used in CloudASM

- bowtie2-2.1.0
- samtools-0.1.19
- bedtools2
- Bismark_v0.19.0
- cutadapt-1.4.2
- TrimGalore-0.4.5
- picard-tools-1.97/picard-tools-1.97
- jre1.7.0_25
- BisSNP-0.82.2

All these packages are included in the Docker-generated image `gcr.io/hackensack-tyco/wgbs-asm`.


## Installation

To run You need to install GCP's Python package called ["dsub"](https://github.com/DataBiosphere/dsub). We recommend using the method where their repository is cloned in combination with the virtual environment.

## How to use the pipeline

1. Prepare the zipped fastq files to be analyzed (see below how)

2. Clone this repository on your computer

3. Define the **ASM variables** and the **GCP variables** in `main.sh`.

4. Launch a virtual environment `source dsub_libs/bin/activate` from dsub's repository.

5. Copy, sequentially, all instructions from main.sh into the terminal, block by block (a "block" is a set of instructions included in between two headlines). 

6. Before moving on to the next instructions block, re-run the jobs that fail if you use preemptive machines (failure rate is about 5-10% when using preemptive machines). See below how to do that.


## Prepare the fastq files to be analyzed


Create a bucket (variable `INPUT_B` in the master script)


All samples you want to analyze need to be in the same bucket (which we call here `gs://SAMPLES`) with one folder per sample. In the bucket, at `gs://SAMPLES/lanes.csv`, upload a CSV file with the correspondance of each zipped fastq file to the lane ID, read, and size in GB of the zipped fastq file. For ENCODE samples, it looks like this:


| sample | bucket_url | lane_id | read_id | file_new_name |
| ------ | ---------- | ------- | ------- | ------------- |
| gm12878 | gs://encode-wgbs/gm12878/ENCFF113KRQ.fastq.gz	| L01 | R2 | gm12878_L01.R2.fastq |
| gm12878 | gs://encode-wgbs/gm12878/ENCFF585BXF.fastq.gz | L02 | R1 | gm12878_L02.R1.fastq |
| gm12878 | gs://encode-wgbs/gm12878/ENCFF798RSS.fastq.gz | L01 | R1 | gm12878_L01.R1.fastq |
| gm12878 | gs://encode-wgbs/gm12878/ENCFF851HAT.fastq.gz | L02 | R2 | gm12878_L02.R2.fastq |





## Re-run failed jobs

```
JOB="variant_call_rerun"
dstat --provider google-v2 --project PROJECT --jobs 'JOB-ID' --users 'USER' --status '*' > JOB.log
cat $JOB.log | grep -v Success | tail -n +3 | awk '{print $2}' > ${JOB}_failed.txt
sed -i '$ d' ${JOB}_failed.txt


head -1 ${JOB}.tsv > ${JOB}_rerun.tsv

while read INDEX ; do
  ROW=$(($INDEX +1))
  sed -n "${ROW}p" ${JOB}.tsv >> ${JOB}_rerun.tsv
done < ${JOB}_failed.txt