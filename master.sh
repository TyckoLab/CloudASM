
#----------------------------------- Library parameters

# These parameters depend on how your libraries were prepared.

## TRIMMING

# Adapters. By default: Illumina adapters in paired-end samples. 
ADAPTER_A="AGATCGGAAGAGCACACGTCTGAAC"
ADAPTER_A2="AGATCGGAAGAGCGTCGTGTAGGGA"


## NET METHYLATION

# Number of base pairs to ignore on each side of each read for net methylation
IGNORE_R1="3" # 5'
THREE_PRIME_IGNORE_R1="3"
IGNORE_R2="2" # 5'
THREE_PRIME_IGNORE_R2="2"

#----------------------------------- ASM parameters

# Reference genome to align the data
GENOME="hg19" # "GRCh38" or "hg19"

# SNP database used to "destroy" CpG sites that overlap with them
SNPS_FOR_CPG="common_snp" # "raw.vcf" or "filtered.vcf" or "common_snp"
SNP_FREQ="0.05" # only used if the option "common_snp" is selected.

# Effect size required at the ASM region level.
ASM_REGION_EFFECT="0.2"

# Minimum CpG coverage required per allele for single CpGs to be considered for CpG ASM, 
CPG_COV="5"

# Minimum number of CpGs we require near a SNP for it to be considered for an ASM region
CPG_PER_ASM_REGION="3"

# In an ASM region, minimum bumber of CpGs with significant ASM in the same direction
CPG_SAME_DIRECTION_ASM="3"

# Number of consecutive CpGs with significant ASM in the same direction (among all well-covered CpGs)
CONSECUTIVE_CPG="2" 

# Minimum reading score of the SNP (in ASCII)
SNP_SCORE="33" # In ASCII, "33" corresponds to a quality score of zero. See https://www.drive5.com/usearch/manual/quality_score.html

# Benjamin-Hochberg threshold
BH_THRESHOLD="0.05"

# p-value cut-off used in all tests for significance
P_VALUE="0.05"


#----------------------------------- GCP parameters 

# GCP global variables
PROJECT_ID="hackensack-tyco"
REGION_ID="us-central1"
ZONE_ID="us-central1-b"

# Name of the Big Query dataset where data will be stored (do not use dashes in the name)
# Use a name to identify your run.
DATASET_ID="cloudasm" 

# Cloud storage variables (use dashes rather than underscores)
# Make sure the buckets exist and are "standard" (as opposed to nearline or coldstorage)
INPUT_B="cloudasm/fastq" # where your zipped fastq are located. Our bucket gs://cloudasm is public with test data in it.
OUTPUT_B="cloudasm" # will be created by the script if it does not exist. You will not be able to write in the public cloudasm public bucket 
REF_DATA_B="wgbs-ref-files" # will be created by the script

# Path of where you downloaded the Github scripts
SCRIPTS="$HOME/GITHUB_REPOS/CloudASM/"

########################## Useful paths (DO NOT MODIFY) ################################

# Path where to store the logs of the jobs
LOG="gs://$OUTPUT_B/logging/${DATASET_ID}"

# Folder to the bisulfite-converted reference genome
REF_GENOME="gs://$REF_DATA_B/$GENOME/ref_genome" 

# Variant database used in SNP calling
ALL_VARIANTS="gs://$REF_DATA_B/$GENOME/variants/*.vcf" 

# Docker image with genomic packages 
DOCKER_GENOMICS="gcr.io/hackensack-tyco/wgbs-asm"

# Light-weight python Docker image with statistical packages.
DOCKER_PYTHON="gcr.io/hackensack-tyco/python"

# Off-the-shelf Docker image for GCP-only jobs
DOCKER_GCP="google/cloud-sdk:255.0.0"

# Create a local folder on the computer 
mkdir -p ${SCRIPTS}/"run_files" 
cd ${SCRIPTS}/"run_files"

# Create dataset in Big Query
bq --location=us mk --dataset ${PROJECT_ID}:${DATASET_ID}


########################## Prepare sample info file ################################

# Refer to the Github README on how to prepare and download the "samples.tsv file below


######################### Prepare generic task files based on the samples #########

# Making sure the sample info file is in the right format and in the "run_files" folder
dos2unix samples.tsv 

#--------- Tasks centered on samples

# List of samples
awk -F "\t" '{if (NR!=1) print $1}' samples.tsv | uniq > sample_id.txt

# Prepare TSV file with just the samples (used for most jobs)
echo -e "--env SAMPLE" > all_samples.tsv

while read SAMPLE ; do
    echo -e "${SAMPLE}" >> all_samples.tsv
done < sample_id.txt

#--------- Tasks centered on chromosomes

# Prepare TSV file per chromosome (used for many jobs)
echo -e "--env SAMPLE\t--env CHR" > all_chr.tsv

# Create a file of job parameters for finding SNPs and their reads.
while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do
      echo -e "${SAMPLE}\t${CHR}" >> all_chr.tsv
  done
done < sample_id.txt


########################## Assemble and prepare the ref genome. Download variants database ################################

# We assemble the ref genome, prepare it to be used by Bismark, and download/unzip the variant database
# This step takes about 6 hours

dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --disk-size 800 \
  --machine-type n1-standard-4 \
  --env GENOME="${GENOME}" \
  --logging $LOG \
  --output-recursive OUTPUT_DIR="gs://${REF_DATA_B}/${GENOME}" \
  --script ${SCRIPTS}/preparation.sh \
  --wait


########################## Unzip, rename, and split fastq files ################################

# Create an TSV file with parameters for the job
echo -e '--input ZIPPED\t--env FASTQ\t--output OUTPUT_FILES' > decompress.tsv

awk -v INPUT_B="${INPUT_B}" \
    -v OUTPUT_B="${OUTPUT_B}" \
    'BEGIN { FS=OFS="\t" } 
    {if (NR!=1) 
        print $2, $5, "gs://"OUTPUT_B"/"$1"/split_fastq/*.fastq" 
     }' \
    samples.tsv >> decompress.tsv 

# Creating ~ 4,000 pairs of 1.2M-row fastq files over 2 hours if zipped fastq file are ~80GB each. 

# Launch job
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --logging $LOG \
  --machine-type n1-standard-2 \
  --disk-size 2000 \
  --preemptible \
  --image $DOCKER_GCP \
  --command 'gunzip ${ZIPPED} && \
             mv ${ZIPPED%.gz} $(dirname "${ZIPPED}")/${FASTQ} && \
             split -l 1200000 \
                --numeric-suffixes --suffix-length=4 \
                --additional-suffix=.fastq \
                $(dirname "${ZIPPED}")/${FASTQ} \
                $(dirname "${OUTPUT_FILES}")/${FASTQ%fastq}' \
  --tasks decompress.tsv \
  --wait

########################## Trim a pair of fastq shards ################################

# Takes about 5 minutes per pair of reads.
# When using preemptive machines, we have experienced a 0.75% failure rate.

# Create an TSV file with parameters for the job
echo -e "--input R1\t--input R2\t--output FOLDER" > trim.tsv

# Prepare inputs and outputs for each sample
while read SAMPLE ; do
  # Get the list of split fastq files
  gsutil ls gs://$OUTPUT_B/$SAMPLE/split_fastq > fastq_shard_${SAMPLE}.txt
  
  # Isolate R1 files
  cat fastq_shard_${SAMPLE}.txt | grep R1 > R1_files_${SAMPLE}.txt && sort R1_files_${SAMPLE}.txt
  # Isolate R2 files
  cat fastq_shard_${SAMPLE}.txt | grep R2 > R2_files_${SAMPLE}.txt && sort R2_files_${SAMPLE}.txt
  # Create a file repeating the output dir for the pair
  NB_PAIRS=$(cat R1_files_${SAMPLE}.txt | wc -l)
  rm -f output_dir_${SAMPLE}.txt && touch output_dir_${SAMPLE}.txt 
  for i in `seq 1 $NB_PAIRS` ; do 
    echo 'gs://'$OUTPUT_B'/'$SAMPLE'/trimmed_fastq/*' >> output_dir_${SAMPLE}.txt
  done
  
  # Add the sample's 3 info (R1, R2, output folder) to the TSV file
  paste -d '\t' R1_files_${SAMPLE}.txt R2_files_${SAMPLE}.txt output_dir_${SAMPLE}.txt >> trim.tsv
done < sample_id.txt

# Print a message in the terminal
echo "There are" $(cat trim.tsv | wc -l) "to be launched"

# Note: you might have to change the adapters and the trimming parameters below based on the library preparation
# By default, we use the Illumina library parameters.

# Submit job. 
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --machine-type n1-standard-2 \
  --preemptible \
  --logging $LOG \
  --env ADAPTER_A="${ADAPTER_A}" \
  --env ADAPTER_A2="${ADAPTER_A2}" \
  --command 'trim_galore \
      -a ${ADAPTER_A} \
      -a2 ${ADAPTER_A2} \
      --quality 30 \
      --length 40 \
      --paired \
      --retain_unpaired \
      --fastqc \
      ${R1} \
      ${R2} \
      --output_dir $(dirname ${FOLDER})' \
  --tasks trim.tsv \
  --wait


########################## Align a pair of fastq shards ################################

# Takes ~10 min per pair of trimmed reads
# About 10% of jobs will fail because GCP will claim back the preemptive machines

# Prepare TSV file
echo -e "--input R1\t--input R2\t--output OUTPUT_DIR" > align.tsv

# Prepare inputs and outputs for each sample
while read SAMPLE ; do
  # Get the list of split fastq files
  gsutil ls gs://$OUTPUT_B/$SAMPLE/trimmed_fastq/*val*.fq > trimmed_fastq_shard_${SAMPLE}.txt
  
  # Isolate R1 files
  cat trimmed_fastq_shard_${SAMPLE}.txt | grep R1 > R1_files_${SAMPLE}.txt && sort R1_files_${SAMPLE}.txt
  # Isolate R2 files
  cat trimmed_fastq_shard_${SAMPLE}.txt | grep R2 > R2_files_${SAMPLE}.txt && sort R2_files_${SAMPLE}.txt
  
  # Create a file repeating the output dir for the pair
  NB_PAIRS=$(cat R1_files_${SAMPLE}.txt | wc -l)
  rm -f output_dir_${SAMPLE}.txt && touch output_dir_${SAMPLE}.txt 
  for i in `seq 1 $NB_PAIRS` ; do 
    echo 'gs://'$OUTPUT_B'/'$SAMPLE'/aligned_per_chard/*' >> output_dir_${SAMPLE}.txt
  done
  
  # Add the sample's 3 info (R1, R2, output folder) to the TSV file
  paste -d '\t' R1_files_${SAMPLE}.txt R2_files_${SAMPLE}.txt output_dir_${SAMPLE}.txt >> align.tsv
done < sample_id.txt

# Print a message in the terminal
echo "There are" $(cat align.tsv | wc -l) "to be launched"

# Submit job (will require about 64,000 CPUs per sample)
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --machine-type n1-standard-16 \
  --preemptible \
  --disk-size 40 \
  --logging $LOG \
  --input-recursive REF_GENOME="${REF_GENOME}" \
  --command 'bismark_nozip \
                -q \
                --bowtie2 \
                ${REF_GENOME} \
                -N 1 \
                -1 ${R1} \
                -2 ${R2} \
                --un \
                --score_min L,0,-0.2 \
                --bam \
                --multicore 3 \
                -o $(dirname ${OUTPUT_DIR})' \
  --tasks align.tsv \
  --wait


########################## Split chard's BAM by chromosome ################################

# Prepare TSV file
echo -e "--input BAM\t--output OUTPUT_DIR" > split_bam.tsv

while read SAMPLE ; do
  gsutil ls gs://$OUTPUT_B/$SAMPLE/aligned_per_chard/*.bam > bam_per_chard_${SAMPLE}.txt
  NB_BAM=$(cat bam_per_chard_${SAMPLE}.txt | wc -l)
  rm -f output_dir_${SAMPLE}.txt && touch output_dir_${SAMPLE}.txt
  for i in `seq 1 $NB_BAM` ; do 
    echo 'gs://'$OUTPUT_B'/'$SAMPLE'/bam_per_chard_and_chr/*' >> output_dir_${SAMPLE}.txt
  done
  paste -d '\t' bam_per_chard_${SAMPLE}.txt output_dir_${SAMPLE}.txt >> split_bam.tsv
done < sample_id.txt

# Submit job
dsub \
  --project $PROJECT_ID \
  --disk-size 30 \
  --preemptible \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --logging $LOG \
  --script ${SCRIPTS}/split_bam.sh \
  --tasks split_bam.tsv \
  --wait


########################## Merge all BAMs by chromosome, clean them ################################

# May take up to 5 hours for the largest chromosomes.

# Prepare TSV file
echo -e "--env SAMPLE\t--env CHR\t--input BAM_FILES\t--output OUTPUT_DIR" > merge_bam.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
  echo -e "${SAMPLE}\t${CHR}\tgs://$OUTPUT_B/$SAMPLE/bam_per_chard_and_chr/*chr${CHR}.bam\tgs://$OUTPUT_B/$SAMPLE/bam_per_chr/*" >> merge_bam.tsv
  done
done < sample_id.txt

# Submit job
dsub \
  --project $PROJECT_ID \
  --machine-type n1-highmem-8 \
  --preemptible \
  --disk-size 120 \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --logging $LOG \
  --script ${SCRIPTS}/merge_bam.sh \
  --tasks merge_bam.tsv \
  --wait


################################# Net methylation call

# May take up to 5 hours for the largest chromosomes

# Prepare TSV file
echo -e "--input BAM\t--output OUTPUT_DIR" > methyl.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
  echo -e "gs://$OUTPUT_B/$SAMPLE/bam_per_chr/${SAMPLE}_chr${CHR}.bam\tgs://$OUTPUT_B/$SAMPLE/net_methyl/*" >> methyl.tsv
  done
done < sample_id.txt

# Note: you may need to customize the number of base pairs you ignore in the net methylation call.

# Launch job
dsub \
  --project $PROJECT_ID \
  --preemptible \
  --machine-type n1-highmem-8 \
  --disk-size 300 \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --logging $LOG \
  --env IGNORE_R1="${IGNORE_R1}" \
  --env THREE_PRIME_IGNORE_R1="${THREE_PRIME_IGNORE_R1}" \
  --env IGNORE_R2="${IGNORE_R2}" \
  --env THREE_PRIME_IGNORE_R2="${THREE_PRIME_IGNORE_R2}" \
  --command 'bismark_methylation_extractor \
                  -p \
                  --no_overlap \
                  --multicore 3 \
                  --merge_non_CpG \
                  --bedGraph \
                  --counts \
                  --report \
                  --buffer_size 48G \
                  --output $(dirname ${OUTPUT_DIR}) \
                  ${BAM} \
                  --ignore ${IGNORE_R1} \
                  --ignore_3prime ${THREE_PRIME_IGNORE_R1} \
                  --ignore_r2 ${IGNORE_R2} \
                  --ignore_3prime_r2 ${THREE_PRIME_IGNORE_R2}' \
  --tasks methyl.tsv \
  --wait



########################## Re-calibrate BAM  ################################

# This step is required by the variant call Bis-SNP

# Takes 5 hours for the largest chromosomes.

# Prepare TSV file
echo -e "--env SAMPLE\t--env CHR\t--input BAM\t--output OUTPUT_DIR" > bam_recalibration.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
  echo -e "$SAMPLE\t$CHR\tgs://$OUTPUT_B/$SAMPLE/bam_per_chr/${SAMPLE}_chr${CHR}.bam\tgs://$OUTPUT_B/$SAMPLE/recal_bam_per_chr/*" >> bam_recalibration.tsv
  done
done < sample_id.txt

# Re-calibrate the BAM files.
dsub \
  --project $PROJECT_ID \
  --machine-type n1-standard-16 \
  --preemptible \
  --disk-size 400 \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --logging $LOG \
  --input REF_GENOME="${REF_GENOME}/*" \
  --input ALL_VARIANTS="${ALL_VARIANTS}" \
  --script ${SCRIPTS}/bam_recalibration.sh \
  --tasks bam_recalibration.tsv \
  --wait


########################## Variant call  ################################

# Takes 6-7 hours for the largest chromosomes. 
# The largest chromosomes (1-5) should probably be run without
# using preemptible

# Prepare TSV file
echo -e "--env SAMPLE\t--env CHR\t--input BAM_BAI\t--output OUTPUT_DIR" > variant_call.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
  echo -e "$SAMPLE\t$CHR\tgs://$OUTPUT_B/$SAMPLE/recal_bam_per_chr/${SAMPLE}_chr${CHR}_recal.ba*\tgs://$OUTPUT_B/$SAMPLE/variants_per_chr/*" >> variant_call.tsv
  done
done < sample_id.txt

# Run tasks
dsub \
  --project $PROJECT_ID \
  --machine-type n1-standard-16 \
  --disk-size 300 \
  --preemptible \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --logging $LOG \
  --input REF_GENOME="${REF_GENOME}/*" \
  --input ALL_VARIANTS="${ALL_VARIANTS}" \
  --script ${SCRIPTS}/variant_call.sh \
  --tasks variant_call.tsv \
  --wait


################################# Export context files in Big Query ##################

# Erase files just in case
while read SAMPLE ; do
  bq rm -f -t ${DATASET_ID}.${SAMPLE}_CpGOB
  bq rm -f -t ${DATASET_ID}.${SAMPLE}_CpGOT
done < sample_id.txt

# Launch job (24 jobs per sample -- 4 samples max if BigQuery limit stays the same)
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_GCP \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --env OUTPUT_B="${OUTPUT_B}" \
  --command 'for STRAND in "OB" "OT" ; do 
               echo "Strand" $STRAND
               bq --location=US load \
               --replace=false \
               --source_format=CSV \
               --skip_leading_rows 1 \
               --field_delimiter "\t" \
               ${DATASET_ID}.${SAMPLE}_CpG${STRAND} \
               gs://${OUTPUT_B}/${SAMPLE}/net_methyl/CpG_${STRAND}_${SAMPLE}_chr${CHR}.txt \
               read_id:STRING,meth_state:STRING,chr:STRING,pos:INTEGER,meth_call:STRING
               done' \
  --tasks all_chr.tsv \
  --name 'export-cpg' \
  --wait

# Append context files and keep CpGs with 10x coverage min 
# Less than 15 minutes per sample.
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_GCP \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --env CPG_COV="${CPG_COV}" \
  --env OUTPUT_B="${OUTPUT_B}" \
  --script ${SCRIPTS}/append_context.sh \
  --tasks all_samples.tsv \
  --wait


########################## Export recal bam to Big Query, clean, and delete from bucket ################################

###### First convert the BAM into SAM

# Prepare TSV file
echo -e "--input BAM\t--output SAM" > bam_to_sam.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
    echo -e "gs://$OUTPUT_B/${SAMPLE}/recal_bam_per_chr/${SAMPLE}_chr${CHR}_recal.bam\tgs://$OUTPUT_B/${SAMPLE}/sam/${SAMPLE}_chr${CHR}_recal.sam" >> bam_to_sam.tsv
  done
done < sample_id.txt


# Create a SAM in the bucket (1h45 for the largest chromosomes)
dsub \
  --project $PROJECT_ID \
  --preemptible \
  --machine-type n1-standard-2 \
  --disk-size 200 \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --logging $LOG \
  --command 'samtools view -o $SAM $BAM' \
  --tasks bam_to_sam.tsv \
  --name 'bam-to-sam' \
  --wait

######## Export all SAM to BigQuery

# Prepare TSV file
echo -e "--env SAMPLE\t--env SAM" > sam_to_bq.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
    echo -e "$SAMPLE\tgs://$OUTPUT_B/${SAMPLE}/sam/${SAMPLE}_chr${CHR}_recal.sam" >> sam_to_bq.tsv
  done
  # Delete existing SAM on big query
  bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_recal_sam_uploaded
done < sample_id.txt

# We append all chromosomes in the same file.
# Takes 2 minutes
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --command 'bq --location=US load \
               --replace=false \
               --source_format=CSV \
               --field_delimiter "\t" \
               --max_bad_records 1000000000 \
               ${DATASET_ID}.${SAMPLE}_recal_sam_uploaded \
               ${SAM} \
               read_id:STRING,flag:INTEGER,chr:STRING,read_start:INTEGER,mapq:INTEGER,cigar:STRING,rnext:STRING,mate_read_start:INTEGER,insert_length:INTEGER,seq:STRING,score:STRING,bismark:STRING,picard_flag:STRING,read_g:STRING,genome_strand:STRING,NM_tag:STRING,meth:STRING,score_before_recal:STRING,read_strand:STRING' \
  --tasks sam_to_bq.tsv \
  --name 'export-sam' \
  --wait

# Clean the SAM on BigQuery
# 1 minute
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --script ${SCRIPTS}/clean_sam.sh \
  --tasks all_samples.tsv \
  --wait

# Delete the SAM files from the bucket (they take a lot of space) 
# and the raw SAM files from Big Query
# Takes 2 minutes
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --command 'gsutil rm ${SAM} && bq rm -f -t ${DATASET_ID}.${SAMPLE}_recal_sam_uploaded' \
  --tasks sam_to_bq.tsv \
  --name 'delete-sam' \
  --wait


########################## Prepare SNP database to destroy CpGs ################################

# This step removes all CpG sites overlapping with a SNP
# If you select common snps, it takes about 2 minutes

dsub \
  --project $PROJECT_ID \
  --ssh \
  --zones $ZONE_ID \
  --machine-type n1-standard-2 \
  --disk-size 40 \
  --image $DOCKER_GENOMICS \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --env OUTPUT_B="${OUTPUT_B}" \
  --env SNPS_FOR_CPG="${SNPS_FOR_CPG}" \
  --env SNP_FREQ="${SNP_FREQ}" \
  --env GENOME="${GENOME}" \
  --script ${SCRIPTS}/snps_for_cpg.sh \
  --tasks all_samples.tsv \
  --wait


########################## Export to BQ and clean the filtered VCF ##################


# We append all chromosomes files in the same file.
# Takes about 4 minutes
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --env OUTPUT_B="${OUTPUT_B}" \
  --command 'bq rm -f -t ${DATASET_ID}.${SAMPLE}_vcf_uploaded \
            && for CHR in `seq 1 22` X Y ; do 
                sleep 1s \
                && bq --location=US load \
                  --replace=false \
                  --source_format=CSV \
                  --field_delimiter "\t" \
                  --skip_leading_rows 117 \
                  ${DATASET_ID}.${SAMPLE}_vcf_uploaded \
                  gs://$OUTPUT_B/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}_filtered.vcf \
                  chr:STRING,pos:STRING,snp_id:STRING,ref:STRING,alt:STRING,qual:FLOAT,filter:STRING,info:STRING,format:STRING,data:STRING
               done' \
  --tasks all_samples.tsv \
  --name 'upld-filt-vcf' \
  --wait

# Clean the VCF -- create temporary tables (one per chr)
# a few minutes
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --script ${SCRIPTS}/clean_vcf.sh \
  --tasks all_samples.tsv \
  --wait


########################## Find the read IDs that overlap the snp ##################


# Delete all files to be replaced
while read SAMPLE ; do
  bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_reads  
done < sample_id.txt

# This is the most intensive step on BigQuery. The current limit is 100 concomitant queries at the same time
# so you cannot run more than 4 samples at the same time without asking GCP to expand your quotas (easy process).
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --script ${SCRIPTS}/reads_overlap_snp.sh \
  --tasks all_chr.tsv \
  --wait

# Merge all chromosome files into a single file per sample 
# and delete the individual chromosome files
# Takes a few minutes
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --command 'bq rm -f -t ${DATASET_ID}.${SAMPLE}_vcf_reads \
            && for CHR in `seq 1 22` X Y ; do
                 bq cp --append_table ${DATASET_ID}.${SAMPLE}_vcf_reads_chr${CHR} ${DATASET_ID}.${SAMPLE}_vcf_reads
                bq rm -f -t ${DATASET_ID}.${SAMPLE}_vcf_reads_chr${CHR}
               done' \
  --tasks all_samples.tsv \
  --name 'app-vcf-reads' \
  --wait

########################## Tag each read with REF or ALT and then
########################## each pair of SNP and CpG     ##################

# We consider the cases where there are 1, 3, and 5 numbers in the CIGAR string
# We leave out the 0.00093% where the CIGAR string has 7 numbers or more
# Note: 0.05% of SNPs are left out when the SNP is at the last position of the read.

# We also remove the reads where the SNP cannot be identified with a minimal score.

# Takes ~2 min
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --env SNP_SCORE="${SNP_SCORE}" \
  --script ${SCRIPTS}/read_genotype.sh \
  --tasks all_samples.tsv \
  --wait

# Tag the pair (CpG, snp) with REF or ALT

# we remove CpG that do not have at least 5x on both ref and alt
# This is very stringent: only ~10% of CpGs are left.

# We also remove CpG where a SNP occurs on the C or G
# This removes 1% of well-covered CpGs.
# Takes about ~ 10 minutes
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env OUTPUT_B="${OUTPUT_B}" \
  --env DATASET_ID="${DATASET_ID}" \
  --env CPG_COV="${CPG_COV}" \
  --script ${SCRIPTS}/cpg_genotype.sh \
  --tasks all_samples.tsv \
  --wait

########################## Calculate ASM at the single CpG level ##################


# Prepare TSV file
echo -e "--input CPG_GENOTYPE\t--output CPG_ASM" > cpg_asm.tsv

while read SAMPLE ; do
    echo -e "gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_cpg_genotype.csv\tgs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_cpg_asm.csv" >> cpg_asm.tsv
done < sample_id.txt

# Takes about one hour (4 CPU-hours)
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --disk-size 50 \
  --machine-type n1-standard-4 \
  --image ${DOCKER_PYTHON} \
  --logging $LOG \
  --script ${SCRIPTS}/asm_single_cpg.py \
  --tasks cpg_asm.tsv \
  --wait


########################## Constitute the ASM regions ##################


# Prepare TSV file
echo -e "--input ASM_REGION\t--output ASM_REGION_PVALUE" > asm_regions.tsv

while read SAMPLE ; do
    echo -e "gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_snp_for_asm_region.json\tgs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_asm_region_pvalue.json" >> asm_regions.tsv
done < sample_id.txt

# Takes three minute
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --env OUTPUT_B="${OUTPUT_B}" \
  --env CPG_PER_ASM_REGION="${CPG_PER_ASM_REGION}" \
  --env P_VALUE="${P_VALUE}" \
  --script ${SCRIPTS}/asm_region.sh \
  --tasks all_samples.tsv \
  --wait


# Compute Wilcoxon's p-value per asm_region between the REF reads and the ALT reads
# Calculate the number of consecutive ASMs in the same direction
# Takes 30 minutes 
dsub \
  --ssh \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --disk-size 100 \
  --machine-type n1-highmem-4 \
  --image ${DOCKER_PYTHON} \
  --logging $LOG \
  --env P_VALUE="${P_VALUE}" \
  --env BH_THRESHOLD="${BH_THRESHOLD}" \
  --script ${SCRIPTS}/asm_region.py \
  --tasks asm_regions.tsv \
  --wait

########################## Provide a final list of ASM regions ##################

# Takes 3 minutes (0.05 CPU-hours)
dsub \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging $LOG \
  --env DATASET_ID="${DATASET_ID}" \
  --env OUTPUT_B="${OUTPUT_B}" \
  --env ASM_REGION_EFFECT="${ASM_REGION_EFFECT}" \
  --env CPG_SAME_DIRECTION_ASM="${CPG_SAME_DIRECTION_ASM}" \
  --env P_VALUE="${P_VALUE}" \
  --env CONSECUTIVE_CPG="${CONSECUTIVE_CPG}" \
  --script ${SCRIPTS}/summary.sh \
  --tasks all_samples.tsv \
  --wait


############################# Delete intermediary files #######################

# Delete intermediary files on BigQuery

while read SAMPLE ; do
  echo "Working on sample " $SAMPLE "..."
  bq rm -f -t ${DATASET_ID}.${SAMPLE}_context_filtered
  bq rm -f -t ${DATASET_ID}.${SAMPLE}_cpg_read_genotype
  bq rm -f -t ${DATASET_ID}.${SAMPLE}_recal_sam_uploaded
  bq rm -f -t ${DATASET_ID}.${SAMPLE}_snps_for_cpg
  bq rm -f -t ${DATASET_ID}.${SAMPLE}_vcf_reads_genotype
  bq rm -f -t ${DATASET_ID}.${SAMPLE}_vcf_reads_for_genotyping
done < sample_id.txt

# Delete splited fastq files to save space on the bucket.
while read SAMPLE ; do
  touch split_deleted_after_alignment.log
  gsutil cp split_deleted_after_alignment.log gs://$OUTPUT_B/$SAMPLE/split_fastq/deleted_after_alignment.log
  gsutil -m rm gs://$OUTPUT_B/$SAMPLE/split_fastq/*.fastq
done < sample_id.txt

# Delete BAM files split per chard and chromosome
while read SAMPLE ; do
  touch bam_deleted.log
  gsutil cp bam_deleted.log gs://$OUTPUT_B/$SAMPLE/bam_per_chard_and_chr/bam_deleted.log
  gsutil -m rm gs://$OUTPUT_B/$SAMPLE/bam_per_chard_and_chr/*.bam
done < sample_id.txt

# Delete non-CpG files context files
while read SAMPLE ; do
  touch delete_noncpg.log
  gsutil cp delete_noncpg.log gs://$OUTPUT_B/$SAMPLE/net_methyl/delete_noncpg.log
  gsutil -m rm gs://$OUTPUT_B/$SAMPLE/net_methyl/Non_CpG*

# Delete reference genome files
gsutil -m rm -r gs://${REF_DATA_B}/*
