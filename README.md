
# CloudASM: a cloud-based, ultra-efficient pipeline for mapping allele-specific DNA methylation

Last updated: October 11, 2019. Please check our preprint on [biorxiv](https://www.biorxiv.org/).

## Table of contents

  - [Overview](#overview)
  - [Biology significance of CloudASM](#biology-significance-of-cloudasm)
  - [Pipeline overview](#pipeline-overview)
  - [Bioinformatics packages used in CloudASM](#bioinformatics-packages-used-in-cloudasm)
  - [Installation](#installation)
  - [How to use the pipeline](#how-to-use-the-pipeline)
  - [Prepare the fastq files to be analyzed](#prepare-the-fastq-files-to-be-analyzed)
  - [Re-run failed jobs](#re-run-failed-jobs)

***********

## Overview

CloudASM is a turnkey pipeline designed to call allele-specific CpG methylation in whole methylomes. It is designed to run on [Google Cloud Platform](https://cloud.google.com/) (GCP). 

Because it leverages the Google's severless data warehouse, all steps specific to ASM-calling (as opposed to alignment, variant calling, and methylation calling) are about ~500x more efficient in terms of CPU-hours when compared to a traditional server-based approach using a cluster.

This pipeline starts from the zipped fastq files hosted on a bucket and outputs a table on BigQuery of all regions in a sample showing allele-specific methylation.

## Biology significance of CloudASM

Our laboratory has a long-standing expertise in studying allele-specific methylation. To make sure our pipeline avoids outputing false positives, we have implemented the following steps:

- we filter out the reads where the confidence in the SNP nucleotide is lower than 30 (variable `SNP_SCORE`)
- we remove CpGs from the context file where the C or G overlap with a SNP found in the unfiltered list of variants identified by BisSNP
- We do not consider for ASM the SNPs that are not within 500 bp of at least a CpG

To catch a "true positive" phenomenon of allele-specific methylation, we use two types of calculation:

1. single CpG level ASM ("CpG ASM") where we estimate if there is ASM on a CpG that is at least 5x covered on both alleles. We use a cut-off p-value of 0.05 on an exact Fisher's test.
2. ASM over a region delimited by CpGs that are found on the same reads as the SNP. The region is delimited by two CpGs showing ASM, must contain at least 3 CpGs (variable `CPG_PER_DMR`), must have at least 2 consecutive CpGs with ASM in the same direction (variable `CONSECUTIVE_CPG`), must have an effect size of at least 20% between the two allele (variable `DMR_EFFECT`). The effect size is the difference in methylation percentage between the two alleles calculated across all CpGs in the DMR. We use a p-value of 0.05 on a Wilcoxon test.


Note:

- choosing `common_snp` is ~10x larger and will destroy a lot of CpGs.

## Pipeline overview

The pipeline follows these steps:

1. Create a bisulfite-converted genome from the reference genome that was chosen.
2. Unzip fastq files and split them into smaller files.
3. Trim each pair of fastq files
4. Align each pair of fastq files
5. Merge BAM files per chromosome and prepare them for variant calling.
6. Net methylation calling
7. Re-calibration of BAM files
8. Variant calling
9. Allele-specific methylation calling

Note that the user can choose the reference genome to align the bisulfite-converted reads. The script automatically fetches the correct database of variants for the reference genome that is selected.

## Bioinformatics packages used in CloudASM

- bowtie2-2.1.0
- samtools-0.1.19
- bedtools2
- Bismark_v0.19.0
- cutadapt-1.4.2
- TrimGalore-0.4.5
- picard-tools-1.97
- jre1.7.0_25
- BisSNP-0.82.2

All these packages are included in the Docker-generated image `gcr.io/hackensack-tyco/wgbs-asm`.

Note #1: we need to use a specific version of JAVA to be able to run BisSNP.

Note #2: Bismark reports sequences it aligns reads in positive strand but it is bisulfite-converted. Bis-SNP reports SNPs in positive strand of the reference genome, which is not bisulfite-converted. Therefore, handling these two datasets together requires additional steps.

## Installation

To run You need to install GCP's Python package called ["dsub"](https://github.com/DataBiosphere/dsub). We recommend using the method where their repository is cloned in combination with the virtual environment.

## How to use the pipeline

1. Prepare the zipped fastq files to be analyzed (see below how)

2. Clone this repository on your computer

3. Define the **ASM variables** and the **GCP variables** in `main.sh`.

4. Launch a virtual environment `source dsub_libs/bin/activate` from dsub's repository.

5. Copy, sequentially, all instructions from main.sh into the terminal, block by block (a "block" is a set of instructions included in between two headlines). 

6. Before moving on to the next instructions block, re-run the jobs that fail if you use preemptible machines (failure rate is about 5-10% when using preemptible machines). See below how to do that.

## Prepare the fastq files to be analyzed

Create a bucket (variable `INPUT_B` in the master script) and create, within this bucket, one folder per sample (the folder should have the name of the sample -- do not use dashes in the name of the sample). In the `INPUT_B` bucket, upload a CSV file that describes the sample's files, using the model below:

| sample | bucket_url | lane_id | read_id | file_new_name |
| ------ | ---------- | ------- | ------- | ------------- |
| gm12878 | gs://encode-wgbs/gm12878/ENCFF113KRQ.fastq.gz	| L01 | R2 | gm12878_L01.R2.fastq |
| gm12878 | gs://encode-wgbs/gm12878/ENCFF585BXF.fastq.gz | L02 | R1 | gm12878_L02.R1.fastq |
| gm12878 | gs://encode-wgbs/gm12878/ENCFF798RSS.fastq.gz | L01 | R1 | gm12878_L01.R1.fastq |
| gm12878 | gs://encode-wgbs/gm12878/ENCFF851HAT.fastq.gz | L02 | R2 | gm12878_L02.R2.fastq |

Do not change the titles of the columns of this CSV file.

## Re-run failed jobs

If you use preemptible machines, they may be taken back by GCP before the job ends. If this is the case, then you need to re-run the tasks that could not be completed before the termination of the preemptible machines executing them. The code below will enable you to create a new TSV of all the tasks that could not complete on time.

```
JOB="TASK"
dstat --provider google-v2 --project PROJECT --jobs 'JOB-ID' --users 'USER' --status '*' > JOB.log
cat $JOB.log | grep -v Success | tail -n +3 | awk '{print $2}' > ${JOB}_failed.txt
sed -i '$ d' ${JOB}_failed.txt

head -1 ${JOB}.tsv > ${JOB}_rerun.tsv

while read INDEX ; do
  ROW=$(($INDEX +1))
  sed -n "${ROW}p" ${JOB}.tsv >> ${JOB}_rerun.tsv
done < ${JOB}_failed.txt
```