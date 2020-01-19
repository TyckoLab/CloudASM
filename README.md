
# CloudASM: a cloud-based, ultra-efficient pipeline for mapping allele-specific DNA methylation

Last updated: October 11, 2019. Please check our preprint on [biorxiv](https://www.biorxiv.org/).

## Table of contents

  - [Overview](#overview)
  - [Biology significance of CloudASM](#biology-significance-of-cloudasm)
  - [Main steps in the CloudASM pipeline](#Main-steps-in-the-CloudASM-pipeline)
  - [Bioinformatics packages used in CloudASM](#bioinformatics-packages-used-in-cloudasm)
  - [Installation](#installation)
  - [How to use the pipeline](#how-to-use-the-pipeline)
  - [Prepare the fastq files to be analyzed](#prepare-the-fastq-files-to-be-analyzed)
  - [Re-run failed jobs](#re-run-failed-jobs)

***********

## Overview

CloudASM is a turnkey pipeline designed to call allele-specific CpG methylation in whole methylomes. It is designed to run on [Google Cloud Platform](https://cloud.google.com/) (GCP). 

This pipeline takes as an input  zipped fastq files and outputs a table of all single nucleotide polymorphisms (SNPs) with allele-specific methylation. Below, we show an example of the output table:

|chr|snp_id|snp_pos|asm_snp|asm_region_inf|asm_region_sup|nb_ref_reads|nb_alt_reads|asm_region_effect|wilcoxon_corr_pvalue|nb_cpg|nb_sig_cpg|nb_pos_sig_cpg|nb_neg_sig_cpg|nb_consec_pos_sig_asm|nb_consec_neg_sig_asm|
|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|
|1|rs1009940|9278142|true|9278133|9278266|15|15|-0.494|0.00322|5|4|1|3|0|2|
|1|rs10127939|161518333|true|161518283|161518415|24|17|-0.297|0.03277|6|3|0|3|0|2|
|1|rs10157041|161397958|true|161397867|161397946|21|35|-0.426|0.00275|12|8|0|8|0|6|

Here are the explanations of the different parameters on this table:

- `chr`: this this the chromosome number where the SNP is located
- `snp_id`: the unique identifier for the SNP that was found to have ASM.
- `snp_pos`: Coordinate of the SNP.
- `asm_region_inf`: Position of the CpG with significant ASM for the SNP `snp_id` and the smallest coordinate.
- `asm_region_sup`: Position of the CpG with significant ASM for the SNP `snp_id` and the smallest coordinate.
- `nb_ref_reads`: Number of genomic segments that cover the REF of the SNP.
- `nb_alt_reads`: Number of genomic segments that cover the ALT of the SNP.
- `asm_region_effect`: The difference between fractional methylation between the ALT and REF genomics segments.
- `wilcoxon_corr_pvalue`: Wilcoxon p-value of the the `asm_region_effect` corrected for multiple testing (across all SNPs).
- `nb_cpg`: Number of CpGs with at least 5x coverage on both the REF and ALT genomic segments.
- `nb_pos_sig_cpg`: Number of CpGs with significant ASM and where fractional methylation between ALT and REF is positive.
- `nb_neg_sig_cpg`: Number of CpGs with significant ASM and where fractional methylation between ALT and REF is negative.
- `nb_consec_pos_sig_asm`: Number of consecutive CpGs with significant ASM and where fractional methylation between ALT and REF is positive.
- `neg_consec_sig_cpg`: Number consecutive of CpGs with significant ASM and where fractional methylation between ALT and REF is negative.

After running the pipeline, the results will be available in a CSV file at `gs://OUTPUT_B/SAMPLE/asm/SAMPLE_asm.csv` where `SAMPLE` is the name of the sample and `OUTPUT_B` is the name of the bucket (one of the variables in the master script).

## Biology significance of CloudASM

Our laboratory has a long-standing expertise in studying allele-specific methylation. To make sure our pipeline avoids outputing false positives, we have implemented the following steps for stringency:

![Definition of an ASM region]([http://url/to/img.png](https://github.com/TyckoLab/CloudASM/blob/master/ASM.png))


- we filter out the reads where the confidence in the SNP nucleotide is lower than 30 (variable `SNP_SCORE`)
- we remove CpGs from the context file where the C or G overlap with a SNP found in the unfiltered list of variants identified by BisSNP
- We do not consider for ASM the SNPs that are not within 500 bp of at least a CpG

To catch a "true positive" phenomenon of allele-specific methylation, we use two types of calculation:

1. single CpG level ASM ("CpG ASM") where we estimate if there is ASM on a CpG that is at least 5x covered on both alleles. We use a cut-off p-value of 0.05 on an exact Fisher's test.
2. ASM over a region delimited by CpGs that are found on the same reads as the SNP. The region is delimited by two CpGs showing ASM, must contain at least 3 CpGs (variable `CPG_PER_asm_region`), must have at least 2 consecutive CpGs with ASM in the same direction (variable `CONSECUTIVE_CPG`), must have an effect size of at least 20% between the two allele (variable `ASM_REGION_EFFECT`). The effect size is the difference in methylation percentage between the two alleles calculated across all CpGs in the asm_region. We use a p-value of 0.05 on a Wilcoxon test.

All these variables can be adjusted by the user.

## Main steps in the CloudASM pipeline

The pipeline follows these steps:

1. Create a bisulfite-converted genome from the reference genome that was chosen
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

All these packages are included in the Docker-generated publicly-available image `gcr.io/hackensack-tyco/wgbs-asm`. Note that we need to use a specific version of JAVA to be able to run BisSNP.

## If you are new to Google Cloud Computing

Google Cloud Computing (GCP) is a service for cloud computing, just like Amazon Web Services and Microsoft Azur Cloud Compute.

To be able use CloudASM, you need to create an account on https://cloud.google.com/ and set up how you want to be billed (e.g. with a credit card). As of January 2020, GCP offers $300 in credits for opening an account -- which is enough to test this pipeline.

Once you do that, you need to create a "project" within GCP and choose which geographical ["region" and "zone"](https://www.google.com/search?q=gcp+regions&rlz=1C5CHFA_enUS809US809&oq=gcp+regions&aqs=chrome..69i57.1415j0j7&sourceid=chrome&ie=UTF-8) you want to request resources from. It is recommended to pick a region and zone near your location. Every job you submit within your project will pull resources from the region and zone you chose.

EXPLAIN HOW THE DATA IS IN "BUCKETS" ETC.

jobs with dsub requires the project name, the zone in which the
virtual machine (VM) is to be launched, input bucket, output
bucket, logging bucket, docker image name, machine configuration and a command/program to execute. Upon receiving this
information, Google cloud executes that command on a docker
image launched on . Upon receiving this
information, Google cloud executes that command on a docker
image launched on a VM with the machine configuration
specified by dsub.


For this reason, you will need to go over to the [Quotas](https://console.cloud.google.com/iam-admin/quotas) and you need to make sure you have access to the following quotas to be able to run 10 WGBS samples at the same time:
- 3,000 `Queries per 100 seconds` in Compute Engine API (zone: "Global")
- 2,000 `Read requests per 100 seconds` in Compute Engine API (zone: "Global")
- 2,000 `Lists requests per 100 seconds` in Compute Engine API (zone: "Global")
- 50,000 `Queries per 100 seconds` for the Genomics API (zone: "Global")
- 2,000 `In-use IP addresses` in Compute Engine API in your zone.
- 200 TB of `Persistent Disk Standard (GB)` in your zone.
- 100,000 CPUs in Compute Engine API in your zone.

All of these values should be by default except the number of CPUs. We need 16 x 


GCP offers a suite of services and CloudASM uses [Compute Engine](https://cloud.google.com/compute/), where virtual machines can be created and used for CloudASM's tasks, [Cloud Storage](https://cloud.google.com/storage/), where data is stored before and after being processed by virtual machines, and [Big Query](https://cloud.google.com/bigquery/), where the aligned reads, variants, and methylation states are analyzed jointly to estimate ASM.

## Installation

To run CloudASM, you need to install GCP's Python package called ["dsub"](https://github.com/DataBiosphere/dsub). We recommend using the method where their repository is cloned in combination with the virtual environment.

## How to use the pipeline

1. Prepare the zipped fastq files to be analyzed (see below how)

2. Clone this repository on your computer

3. Define the **ASM variables** and the **GCP variables** in `master.sh`.

4. Launch a virtual environment `source dsub_libs/bin/activate` from dsub's repository.

5. Copy, sequentially, all instructions from `master.sh` into the terminal, block by block (a "block" is a set of instructions included in between two headlines). 

6. Before moving on to the next instructions block, re-run the jobs that fail if you use preemptible machines (failure rate is about 5-10% when using preemptible machines). See below how to do that.

## Test the pipeline

https://console.cloud.google.com/storage/browser/cloudasm

## Prepare the fastq files to be analyzed

https://docs.google.com/spreadsheets/d/1WFpR_uM9BdBAdoCcIoVfM7rjuJF-khlvJTKHHd7bkEQ/edit#gid=0

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