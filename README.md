
# CloudASM: an ultra-efficient cloud-based pipeline for mapping allele-specific DNA methylation

Last updated: July 31, 2020. Please check our preprint on [biorxiv](https://www.biorxiv.org/content/10.1101/2020.01.28.887430v1) and our published paper on [Bioinformatics](https://doi.org/10.1093/bioinformatics/btaa149).

## Table of contents

  - [Overview](#overview)
  - [Biology significance of CloudASM](#biology-significance-of-cloudasm)
  - [Main steps in the CloudASM pipeline](#Main-steps-in-the-CloudASM-pipeline)
  - [Bioinformatics packages used in CloudASM](#bioinformatics-packages-used-in-cloudasm)
  - [If you are new to Google Cloud Computing](#if-you-are-new-to-google-cloud-computing)
  - [Installation](#installation)
  - [How to use the pipeline](#how-to-use-the-pipeline)
  - [Test the pipeline](#test-the-pipeline)
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
- `snp_id`: the unique identifier for the SNP that was evaluated for ASM.
- `snp_pos`: Coordinate of the SNP.
- `asm_snp`: whether the SNP has ASM or not.
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

See below how these look like on an example.

## Biology significance of CloudASM

![Definition of an ASM region](https://github.com/TyckoLab/CloudASM/blob/master/ASM.png)

Our laboratory has a long-standing expertise in studying allele-specific methylation. To make sure our pipeline avoids outputing false positives, we have implemented the following steps for stringency:

- we filter out the reads where the confidence in the SNP nucleotide is lower than 30 (variable `SNP_SCORE`)
- we remove CpGs from the context file where the C or G overlap with a SNP found in the unfiltered list of variants identified by BisSNP
- We do not consider for ASM the SNPs that are not within 500 bp of at least a CpG

As described in the figure above, to catch a "true positive" phenomenon of allele-specific methylation, we use two types of calculation:

1. single CpG level ASM ("CpG ASM") where we estimate if there is ASM on a CpG that is at least 5x covered on both alleles. We use a cut-off p-value of 0.05 on an exact Fisher's test.
2. ASM over a region delimited by CpGs that are found on the same reads as the SNP. The region is delimited by two CpGs showing ASM, must contain at least 3 CpGs (variable `CPG_PER_asm_region`), must have at least 2 consecutive CpGs with ASM in the same direction (variable `CONSECUTIVE_CPG`), must have an effect size of at least 20% between the two allele (variable `ASM_REGION_EFFECT`). The effect size is the difference in methylation percentage between the two alleles calculated across all CpGs in the asm_region. We use a p-value of 0.05 on a Wilcoxon test.

All these variables can be adjusted by the user.

## Main steps in the CloudASM pipeline

The pipeline follows these steps:

1. Create a bisulfite-converted genome from the reference genome that was chosen
2. Unzip fastq files and split them into smaller files ("chards").
3. Trim each pair of fastq files
4. Align each pair of fastq files
5. Merge BAM files per chromosome and prepare them for variant and net methylation calling.
6. Net methylation calling
7. Re-calibrate of BAM files
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

GCP offers a suite of services and CloudASM uses [Compute Engine](https://cloud.google.com/compute/), where virtual machines can be created and used for CloudASM's tasks, [Cloud Storage](https://cloud.google.com/storage/), where data is stored before and after being processed by virtual machines, and [Big Query](https://cloud.google.com/bigquery/), where the aligned reads, variants, and methylation states are analyzed jointly to estimate ASM.

Once you open an account on GCP, you need to create a "project" within GCP and choose which geographical ["region" and "zone"](https://www.google.com/search?q=gcp+regions&rlz=1C5CHFA_enUS809US809&oq=gcp+regions&aqs=chrome..69i57.1415j0j7&sourceid=chrome&ie=UTF-8) you want to request resources from. It is recommended to pick a region and zone near your location. Every job you submit within your project will pull resources from the region and zone you chose.

Very basically, data is stored in "buckets" in the module `Cloud Storage`. When CloudASM executes a batch of jobs, virtual machines (also called "instances") in the module `Compute Engine` are created to execute the job. They download the data from the bucket, obtain the script from CloudASM, execute the script, and upload the output of the script back in the bucket (or BigQuery). Depending on the job, CloudASM requests instances with 2-16 CPUs, adds larger disks to these instances (30-400 GB), and executes the commands on one of the Docker images we built or that were already available.

When you start running CloudASM on more than one sample, pipeline manager dsub will start queuing jobs in the same batch if you do not have enough resoures for your project. When you start running CloudASM, you may want to go to the [Quotas](https://console.cloud.google.com/iam-admin/quotas) and monitor which resources you need to increase.

## Installation

To run CloudASM, you need to install GCP's pipeline manager called ["dsub"](https://github.com/DataBiosphere/dsub). We recommend using the method where their repository is cloned in combination with the virtual environment. You can also [install dsub via conda](https://anaconda.org/conda-forge/dsub) but this was not created by the dsub creators so it may not have the latest version of dsub.

The current version of CloudASM was validated with dsub 0.3.7.

## How to use the pipeline

1. Prepare the zipped fastq files to be analyzed (see below how)

2. Clone this repository on your computer

3. Customize the **Library parameters**, the **ASM parameters** and the **GCP parameters** in `master.sh`.

4. Launch a virtual environment by typing `python3 -m venv dsub_libs` and then `source dsub_libs/bin/activate` (does not matter in which directory you are located).

5. Copy, sequentially, all instructions from `master.sh` into the terminal, block by block (a "block" is a set of instructions included in between two headlines). 

6. Before moving on to the next instructions block, re-run the jobs that fail if you use preemptible machines (failure rate is about 5-10% when using preemptible machines). See below how to do that.

## Test the pipeline

We prepared a small dataset (~24MB of zipped fastq files) for anyone to test CloudASM without incurring large costs. The dataset was prepared by recombining bisulfite-converted reads overlapping the positions 9,000,000 and 10,000,000 on chromosome 1, using the lymphoblastoid cell line GM12878 (identifier: ENCSR890UQO) made publicly available by the ENCODE consortium.

The zipped fastq files are freely accessible on [our GCP’s storage space](https://console.cloud.google.com/storage/browser/cloudasm/fastq/gm12878/). Using this dataset, CloudASM assessed 456 SNPs and found 13 ASM regions. All the data generated by CloudASM for this dataset is stored [here](https://console.cloud.google.com/storage/browser/cloudasm). The final table of all variants and their assessment for ASM can be [downloaded here](https://storage.cloud.google.com/cloudasm/gm12878/asm/gm12878_asm.csv).

To test the pipeline, you will have to change all GCP parameters except the variable `INPUT_B`. As you will when running the pipeline, for each sample, CloudASM creates the following folders:

- `split_fastq`: fastq files unzipped and split into 120,000 rows chards
- `trimmed_fastq`: trimmed pairs of fastq chards from `split_fastq`
- `aligned_per_chard`: BAM files created from aligning trimmed chards of fastq files from `trimmed_fastq`
- `bam_per_chard_and_chr`: The BAM files from `aligned_per_chard` are split across chromosomes
- `bam_per_chr`: The BAM files from `bam_per_chard_and_chr` are merged across chromosomes
- `net_methyl`: context files created from calling net methylation on the BAM files located in `bam_per_chr`.
- `recal_bam_per_chr`: Recalibrated BAM files from `bam_per_chr` 
- `variants_per_chr`: VCF files for each chromosomes created from the BAM files in `recal_bam_per_chr`.
- `asm`: Tables of CpG-level and region-level of ASM.
- `bedgraph`: bedgraph files of coverage and methylation percentage across the whole sample
- `sam`: SAM files created from BAM files in `recal_bam_per_chr`.
- `snps_for_cpg`: database of SNPs used to remove CpG sites potentially overlapping with a variant.

## Prepare the fastq files to be analyzed

When naming your samples, do not use dashes. Only use underscores.

Prepare a sample info file using [this model](https://docs.google.com/spreadsheets/d/1WFpR_uM9BdBAdoCcIoVfM7rjuJF-khlvJTKHHd7bkEQ/edit#gid=0). Download the file as TSV file into your `run_files` directory located where you cloned the CloudASM repository (`run_files` is automatically created by the `master` script). The sample info file looks like the table below. The first column is the name of the sample, the 2nd column is the list of all the zipped fastq files, the 3rd column is the lane of the zipped fastq files, the 4th column is the read number of the zipped fastq file. The 4th column is created automatically from column 1, 3, and 4. 

| sample | bucket_url | lane_id | read_id | file_new_name |
| ------ | ---------- | ------- | ------- | ------------- |
| gm12878 | gs://encode-wgbs/gm12878/ENCFF113KRQ.fastq.gz	| L01 | R2 | gm12878_L01.R2.fastq |
| gm12878 | gs://encode-wgbs/gm12878/ENCFF585BXF.fastq.gz | L02 | R1 | gm12878_L02.R1.fastq |
| gm12878 | gs://encode-wgbs/gm12878/ENCFF798RSS.fastq.gz | L01 | R1 | gm12878_L01.R1.fastq |
| gm12878 | gs://encode-wgbs/gm12878/ENCFF851HAT.fastq.gz | L02 | R2 | gm12878_L02.R2.fastq |
