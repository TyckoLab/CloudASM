
########################## ASM Variables ################################

# Reference genome
GENOME="hg19" # or "GRCh38"

# Effect size required at the DMR level for an ASM.
DMR_EFFECT="0.2"

# Minimum CpG coverage required per allele for single CpGs to be considered for CpG ASM, in a DMR, or "near" a SNP
CPG_COV="5"

# Number of CpGs we require near a SNP for it to be tested with ASM DMR 
# In a DMR, it is also the number of CpGs with significant ASM in the same direction
CPG_PER_DMR="3"

# Number of consecutive CpGs with significant ASM in the same direction (among all well-covered CpGs)
CONSECUTIVE_CPG="2" 

# Minimum reading score of the SNP
SNP_SCORE="63" # In ASCII, 63 correponds to a quality score of 30. See this table: https://www.drive5.com/usearch/manual/quality_score.html

########################## GCP variables ################################

# GCP global variables
PROJECT_ID="hackensack-tyco"
REGION_ID="us-central1"
ZONE_ID="us-central1-b"

# Big Query variables
DATASET_ID="wgbs_encode" 

# Cloud storage variables
INPUT_B="encode-wgbs"
OUTPUT_B="em-encode-paper" # will be created by the script
REF_DATA_B="wgbs-ref-files" 
REF_GENOME="gs://$REF_DATA_B/$GENOME/ref_genome"
ALL_VARIANTS="gs://$REF_DATA_B/$GENOME/variants/*.vcf"

# Path of where you downloaded the Github scripts
SCRIPTS="$HOME/GITHUB_REPOS/wgbs-asm/"

########################## Docker Images used by the pipeline ################################

# Genomic packages 
DOCKER_GENOMICS="gcr.io/hackensack-tyco/wgbs-asm"

# Light-weight python image with statistical packages.
DOCKER_PYTHON="gcr.io/hackensack-tyco/python"

# Off-the-shelf GCP image for GCP-only jobs
DOCKER_GCP="google/cloud-sdk:255.0.0"


########################## Download sample info file ################################

# Refer to the Github on how to prepare the sample info file

# Create a local folder on the computer 
mkdir -p $HOME/"wgbs" 
cd $HOME/"wgbs" 

# Download the metadata to the local folder
gsutil cp gs://$INPUT_B/samples.tsv $HOME/"wgbs"
dos2unix samples.tsv 

# List of samples
awk -F "\t" \
    '{if (NR!=1) \
    print $1}' samples.tsv | uniq > sample_id.txt

echo "There are" $(cat sample_id.txt | wc -l) "samples to be analyzed"

# Prepare TSV file with just the samples (used for most jobs)
echo -e "--env SAMPLE" > all_samples.tsv

while read SAMPLE ; do
    echo -e "${SAMPLE}" >> all_samples.tsv
done < sample_id.txt

# Prepare TSV file per chromosome (used for many jobs)
echo -e "--env SAMPLE\t--env CHR" > all_chr.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do
    echo -e "${SAMPLE}\t${CHR}" >> all_chr.tsv
  done
done < sample_id.txt


########################## Create buckets, datasets, and sample info file ################################

# Create a dataset on BigQuery for the samples to be analyzed for ASM
#(Note: very few regions are available for Big Query datasets)
bq --location=us mk --dataset ${PROJECT_ID}:${DATASET_ID}

# Create buckets for the analysis and for the ref genome / variant database
gsutil mb -c standard -l $REGION_ID gs://${OUTPUT_B} 
gsutil mb -c standard -l $REGION_ID gs://${REF_DATA_B}

########################## Assemble and prepare the ref genome. Download variants database ################################

# We assemble the ref genome, prepare it to be used by Bismark, and download/unzip the variant database
# This step takes about 6 hours

# Do it only once! 

dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --disk-size 800 \
  --machine-type n1-standard-4 \
  --env GENOME="${GENOME}" \
  --logging gs://$OUTPUT_B/logging/ \
  --output-recursive OUTPUT_DIR="gs://${REF_DATA_B}/${GENOME}" \
  --script ${SCRIPTS}/preparation.sh \
  --wait


########################## Unzip, rename, and split fastq files ################################

# Takes ~2 hours per 80G of zipped fastq file.

# Create an TSV file with parameters for the job
echo -e '--input ZIPPED\t--env FASTQ\t--output OUTPUT_FILES' > decompress.tsv

awk -v INPUT_B="${INPUT_B}" \
    -v OUTPUT_B="${OUTPUT_B}" \
    'BEGIN { FS=OFS="\t" } 
    {if (NR!=1) 
        print $2, $5, "gs://"OUTPUT_B"/"$1"/split_fastq/*.fastq" 
     }' \
    samples.tsv >> decompress.tsv 

# Creating ~ 4,000 pairs of 1.2M-row fastq files if the zipped fastq file is ~80GB.

# Launch job
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --logging gs://$OUTPUT_B/logging/ \
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
rm -f trim.tsv && touch trim.tsv

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

# Add headers to the file
sed -i '1i --input R1\t--input R2\t--output FOLDER' trim.tsv

# Print a message in the terminal
echo "There are" $(cat trim.tsv | wc -l) "to be launched"

# Submit job. 
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --machine-type n1-standard-2 \
  --preemptible \
  --logging gs://$OUTPUT_B/logging/ \
  --command 'trim_galore \
      -a AGATCGGAAGAGCACACGTCTGAAC \
      -a2 AGATCGGAAGAGCGTCGTGTAGGGA \
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
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --machine-type n1-standard-16 \
  --preemptible \
  --disk-size 40 \
  --logging gs://$OUTPUT_B/logging/ \
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
  --provider google-v2 \
  --project $PROJECT_ID \
  --disk-size 30 \
  --preemptible \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --logging gs://$OUTPUT_B/logging/ \
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
  --provider google-v2 \
  --project $PROJECT_ID \
  --machine-type n1-highmem-8 \
  --preemptible \
  --disk-size 120 \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --logging gs://$OUTPUT_B/logging/ \
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

# Launch job
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --preemptible \
  --machine-type n1-highmem-8 \
  --disk-size 300 \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --logging gs://$OUTPUT_B/logging/ \
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
                  --ignore 3 \
                  --ignore_3prime 3 \
                  --ignore_r2 2 \
                  --ignore_3prime_r2 2' \
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
  --provider google-v2 \
  --project $PROJECT_ID \
  --machine-type n1-standard-16 \
  --preemptible \
  --disk-size 400 \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --logging gs://$OUTPUT_B/logging/ \
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
  for CHR in `seq 12 22` X Y ; do 
  echo -e "$SAMPLE\t$CHR\tgs://$OUTPUT_B/$SAMPLE/recal_bam_per_chr/${SAMPLE}_chr${CHR}_recal.ba*\tgs://$OUTPUT_B/$SAMPLE/variants_per_chr/*" >> variant_call.tsv
  done
done < sample_id.txt

# Run tasks
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --machine-type n1-standard-16 \
  --disk-size 400 \
  --preemptible \
  --zones $ZONE_ID \
  --image $DOCKER_GENOMICS \
  --logging gs://$OUTPUT_B/logging/ \
  --input REF_GENOME="${REF_GENOME}/*" \
  --input ALL_VARIANTS="${ALL_VARIANTS}" \
  --script ${SCRIPTS}/variant_call.sh \
  --tasks variant.tsv \
  --wait


################################# Export context files in Big Query ##################

# Prepare TSV file
echo -e "--env SAMPLE\t--env STRAND\t--env CONTEXT" > context_to_bq.tsv

while read SAMPLE ; do
  for STRAND in OB OT ; do
    for CHR in `seq 1 22` X Y ; do     
        echo -e "${SAMPLE}\t${STRAND}\tgs://$OUTPUT_B/${SAMPLE}/net_methyl/CpG_${STRAND}_${SAMPLE}_chr${CHR}.txt" >> context_to_bq.tsv
    done
  
  # Delete existing context file on big query
  bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_CpG${STRAND}
  done
done < sample_id.txt

# We append all chromosomes in the same file. 
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_GCP \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --command 'bq --location=US load \
               --replace=false \
               --source_format=CSV \
               --skip_leading_rows 1 \
               --field_delimiter "\t" \
               ${DATASET_ID}.${SAMPLE}_CpG${STRAND} \
               ${CONTEXT} \
               read_id:STRING,meth_state:STRING,chr:STRING,pos:INTEGER,meth_call:STRING' \
  --tasks context_to_bq.tsv \
  --name 'export-cpg' \
  --wait


################################# Append context files and keep CpGs with 10x coverage min ########

# Also this script exports bedgraph files of methylation and coverage to look into IGV for instance

# Prepare TSV file
echo -e "--env SAMPLE" > append_context.tsv

while read SAMPLE ; do
  echo -e "$SAMPLE" >> append_context.tsv
done < sample_id.txt

dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_GCP \
  --logging gs://$OUTPUT_B/logging/ \
  --env PROJECT_ID="${PROJECT_ID}" \
  --env DATASET_ID="${DATASET_ID}" \
  --env CPG_COV="${CPG_COV}" \
  --env OUTPUT_B="$OUTPUT_B" \
  --script ${SCRIPTS}/append_context.sh \
  --tasks append_context.tsv \
  --wait



########################## Export recal bam to Big Query, clean, and delete from bucket ################################

# The SAM file was created by the variant_call script
# we do it before cleaning the variant call because it removes the SAM from the bucket, 
# which takes a lot of splace.

# Prepare TSV file
echo -e "--env SAMPLE\t--env SAM" > sam_to_bq.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
    echo -e "$SAMPLE\tgs://$OUTPUT_B/${SAMPLE}/variants_per_chr/${SAMPLE}_chr${CHR}_recal.sam" >> sam_to_bq.tsv
  done
  
  # Delete existing SAM on big query
  bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_recal_sam_uploaded

done < sample_id.txt

# We append all chromosomes in the same file.
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --command 'bq --location=US load \
               --replace=false \
               --source_format=CSV \
               --field_delimiter "\t" \
               --max_bad_records 10 \
               ${DATASET_ID}.${SAMPLE}_recal_sam_uploaded \
               ${SAM} \
               read_id:STRING,flag:INTEGER,chr:STRING,read_start:INTEGER,mapq:INTEGER,cigar:STRING,rnext:STRING,mate_read_start:INTEGER,insert_length:INTEGER,seq:STRING,score:STRING,bismark:STRING,picard_flag:STRING,read_g:STRING,genome_strand:STRING,NM_tag:STRING,meth:STRING,score_before_recal:STRING,read_strand:STRING' \
  --tasks sam_to_bq.tsv \
  --name 'export-sam' \
  --wait

# Clean the SAM on BigQuery
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --env PROJECT_ID="${PROJECT_ID}" \
  --script ${SCRIPTS}/clean_sam.sh \
  --tasks all_samples.tsv \
  --wait

# Delete the SAM files from the bucket (they take a lot of space) 
# and the raw SAM files from Big Query
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --preemptible \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --command 'gsutil rm ${SAM} && bq rm -f -t ${DATASET_ID}.${SAMPLE}_recal_sam_uploaded' \
  --tasks sam_to_bq.tsv \
  --wait


########################## Export raw VCF to Big Query (one per sample) ################################

# The raw VCFs are used to remove CpGs that overlap with a SNP in the raw VCF, confirmed by the variant database.

# Prepare TSV file
echo -e "--env SAMPLE\t--env VCF" > raw_vcf_to_bq.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
    echo -e "$SAMPLE\tgs://$OUTPUT_B/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}_raw.vcf" >> raw_vcf_to_bq.tsv
  done
  
  # Delete existing SAM on big query
  bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_raw_uploaded
done < sample_id.txt

# We append all chromosomes files in the same file.
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --command 'bq --location=US load \
               --replace=false \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 116 \
               ${DATASET_ID}.${SAMPLE}_vcf_raw_uploaded \
               ${VCF} \
               chr:STRING,pos:STRING,snp_id:STRING,ref:STRING,alt:STRING,qual:FLOAT,filter:STRING,info:STRING,format:STRING,data:STRING' \
  --tasks raw_vcf_to_bq.tsv \
  --wait

# Clean the raw VCF that was uploaded and delete the "uploaded" table 
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --script ${SCRIPTS}/clean_raw_vcf.sh \
  --tasks all_samples.tsv \
  --wait



########################## Export to BQ and clean the filtered VCF ##################

# Filter out the SNPs that are not within 500bp of a CpG that is at least 10x covered.
# This removes about 5% of SNPs. Takes ~30min for the small chr.
# Finds out where the SNP can be found with no ambiguity (CT, GA, or both)

# Prepare TSV file
echo -e "--env SAMPLE\t--env VCF" > vcf_to_bq.tsv

while read SAMPLE ; do
  
  for CHR in `seq 1 22` X Y ; do 
    echo -e "$SAMPLE\tgs://$OUTPUT_B/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}_filtered.vcf" >> vcf_to_bq.tsv
  done
  
  # Delete existing VCF files on Big Query
  bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_uploaded
  bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf

done < sample_id.txt

# We append all chromosomes files in the same file.
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --command 'bq --location=US load \
               --replace=false \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 117 \
               ${DATASET_ID}.${SAMPLE}_vcf_filtered_uploaded \
               ${VCF} \
               chr:STRING,pos:STRING,snp_id:STRING,ref:STRING,alt:STRING,qual:FLOAT,filter:STRING,info:STRING,format:STRING,data:STRING' \
  --tasks vcf_to_bq.tsv \
  --wait

# Clean the VCF -- create temporary tables (one per chr)
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --env PROJECT_ID="${PROJECT_ID}" \
  --script ${SCRIPTS}/clean_vcf.sh \
  --tasks all_chr.tsv \
  --wait

# Append all temporary VCFs together into one single clean VCF
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --env PROJECT_ID="${PROJECT_ID}" \
  --command 'bq cp \
               --append_table \
               ${DATASET_ID}.${SAMPLE}_vcf_chr${CHR}_tmp \
               ${DATASET_ID}.${SAMPLE}_vcf \
             && bq rm -f -t ${DATASET_ID}.${SAMPLE}_vcf_chr${CHR}_tmp' \
  --tasks all_chr.tsv \
  --wait

# Delete all raw VCF files from BigQuery
awk 'BEGIN { FS=OFS="\t" } {if (NR!=1) print $0}' all_samples.tsv > just_samples.txt
while read SAMPLE ; do 
  bq rm -f -t ${DATASET_ID}.${SAMPLE}_vcf_raw_uploaded
done < just_samples.txt

########################## Find the read IDs that overlap the snp ##################


# Delete all files to be replaced
while read SAMPLE ; do
  bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_reads  

done < sample_id.txt

# Create one file per chromosome
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --env PROJECT_ID="${PROJECT_ID}" \
  --script ${SCRIPTS}/reads_overlap_snp.sh \
  --tasks all_chr.tsv \
  --wait

# Merge all chromosome files into a single file per sample 
# and delete the individual chromosome files
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --env PROJECT_ID="${PROJECT_ID}" \
  --command 'bq cp \
               --append_table \
               ${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_${CHR} \
               ${DATASET_ID}.${SAMPLE}_vcf_reads \
             && bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_${CHR}' \
  --tasks all_chr.tsv \
  --wait



########################## Tag each read with REF or ALT and then
########################## each pair of SNP and CpG     ##################

# We consider the cases where there are 1, 3, and 5 numbers in the CIGAR string
# We leave out the 0.00093% where the CIGAR string has 7 numbers or more
# Note: 0.05% of SNPs are left out when the SNP is at the last position of the read.

# We also remove the reads where the score of the nucleotide with the SNP is below 30
# This removes ~ 7% of the reads.

# Tag the read
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env OUTPUT_B="${OUTPUT_B}" \
  --env DATASET_ID="${DATASET_ID}" \
  --env PROJECT_ID="${PROJECT_ID}" \
  --env SNP_SCORE="${SNP_SCORE}"
  --script ${SCRIPTS}/read_genotype.sh \
  --tasks all_samples.tsv \
  --wait

# Tag the pair (CpG, snp) with REF or ALT

# we remove CpG that do not have at least 5x on both ref and alt
# This is very stringent: only ~10% of CpGs are left.

# We also remove CpG where a SNP occurs on the C or G
# This removes 1% of well-covered CpGs.
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env OUTPUT_B="${OUTPUT_B}" \
  --env DATASET_ID="${DATASET_ID}" \
  --env PROJECT_ID="${PROJECT_ID}" \
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

# Using python environment.
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --disk-size 50 \
  --machine-type n1-standard-4 \
  --image ${DOCKER_PYTHON} \
  --logging gs://$OUTPUT_B/logging/ \
  --script ${SCRIPTS}/asm_single_cpg.py \
  --tasks cpg_asm.tsv \
  --wait


########################## Constitute the DMRs ##################


# Prepare TSV file
echo -e "--input DMR\t--output DMR_PVALUE" > dmr.tsv

while read SAMPLE ; do
    echo -e "gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_snp_for_dmr.json\tgs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_dmr_pvalue.json" >> dmr.tsv
done < sample_id.txt

# Requesting 
# 1/ at least NB_CPG_SIG CpGs with significative ASM and in the same direction per DMR 
# 2/ a DMR_EFFECT difference between the REF and ALT (computed across reads) across CpGs that are located in between the two extreme significant CpGs
# Removes ~50% of SNPs.
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --env OUTPUT_B="${OUTPUT_B}" \
  --env CPG_PER_DMR="${NB_CPG_SIG}" \
  --script ${SCRIPTS}/dmr.sh \
  --tasks all_samples.tsv \
  --wait


# Compute Wilcoxon's p-value per DMR between the REF reads and the ALT reads
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --disk-size 30 \
  --machine-type n1-standard-4 \
  --image ${DOCKER_PYTHON} \
  --logging gs://$OUTPUT_B/logging/ \
  --script ${SCRIPTS}/dmr.py \
  --tasks dmr.tsv \
  --wait


########################## Provide a final list of DMRs ##################

# Delete master file of all ASM across all samples
while read SAMPLE ; do
  bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.asm_snp

done < sample_id.txt

dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image ${DOCKER_GCP} \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --env OUTPUT_B="${OUTPUT_B}" \
  --env DMR_EFFECT="${DMR_EFFECT}" \
  --script ${SCRIPTS}/summary.sh \
  --tasks all_samples.tsv \
  --wait


#################################################################




# Delete split fastq files to save space on the bucket.
while read SAMPLE ; do
  touch split_deleted_after_alignment.log
  gsutil cp split_deleted_after_alignment.log gs://$OUTPUT_B/$SAMPLE/split_fastq/deleted_after_alignment.log
  gsutil rm gs://$OUTPUT_B/$SAMPLE/split_fastq/*.fastq
done < sample_id.txt
