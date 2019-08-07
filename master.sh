# Pre-requisities


##### Data
# Have a bucket already created with one folder per sample. 
# Each folder should include the zipped fastq files
# In the bucket, there should be a tsv file "samples.tsv" with each file name, its lane number, its read number (R1 or R2), and the size of the zipped fastq in GB

##### Authentification
# Be authentificated by gcloud (type "gcloud config list" to confirm)


########################## Copy and paste within the bash ################################

# GCP variables
PROJECT_ID="hackensack-tyco"
REGION_ID="us-central1"
ZONE_ID="us-central1-b"
DOCKER_IMAGE="gcr.io/hackensack-tyco/wgbs-asm"

# Names of the buckets (must be unique)
INPUT_B="encode-wgbs" # should already exist with the zipped fastq files inside
OUTPUT_B="em-encode-deux" # will be created by the script
REF_DATA_B="wgbs-ref-files" # See documentation for what it needs to contain

# Create a bucket with the analysis
gsutil mb -c regional -l $REGION_ID gs://$OUTPUT_B 

# Create a local directory to download files
WD="$HOME/wgbs" && mkdir -p $WD && cd $WD

# Download the meta information about the samples and files to be analyzed.
gsutil cp gs://$INPUT_B/samples.tsv $WD
dos2unix samples.tsv 

# Check the number of samples to be analyzed.
echo "There are" $(($(awk -F '\t' '{print $1}' samples.tsv | uniq -c | wc -l) - 1)) "samples to be analyzed"

########################## Unzip, rename, and split fastq files ################################

# Create an TSV file with parameters for the job
awk -v INPUT_B="${INPUT_B}" \
    -v OUTPUT_B="${OUTPUT_B}" \
    'BEGIN { FS=OFS="\t" } \
    {if (NR!=1) \
        print "gs://"INPUT_B"/"$1"/"$2, \
              $5".fastq", \
              "gs://"OUTPUT_B"/"$1"/split_fastq/*.fastq" \
     }' \
    samples.tsv > decompress.tsv 

# Add headers to the file
sed -i '1i --input ZIPPED\t--env FASTQ\t--output OUTPUT_FILES' decompress.tsv

# Launch job
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --logging gs://$OUTPUT_B/logging/ \
  --disk-size 10 \
  --preemptible \
  --image $DOCKER_IMAGE \
  --command 'gunzip ${ZIPPED} && \
             mv ${ZIPPED%.gz} $(dirname "${ZIPPED}")/${FASTQ} && \
             split -l 120000 \
                --numeric-suffixes --suffix-length=5 \
                --additional-suffix=.fastq \
                $(dirname "${ZIPPED}")/${FASTQ} \
                $(dirname "${OUTPUT_FILES}")/${FASTQ%fastq}' \
  --tasks decompress.tsv \
  --wait

########################## Trim a pair of fastq shards ################################

# Prepare a TSV file with parameters for the job

# Create an TSV file with parameters for the job
awk -v INPUT_B="${INPUT_B}" \
    -v OUTPUT_B="${OUTPUT_B}" \
    'BEGIN { FS=OFS="\t" } \
    {if (NR!=1) \
        "gs://$OUTPUT_B/A549-extract/split_fastq/A549-extract_L01.R1.001.fastq"
        print "gs://"INPUT_B"/"$1"/"$2, \
              $5".fastq", \
              "gs://"OUTPUT_B"/"$1"/split_fastq/*.fastq" \
     }' \
    samples.tsv > decompress.tsv 

# Add headers to the file
sed -i '1i --input ZIPPED\t--env FASTQ\t--output OUTPUT_FILES' decompress.tsv


# Column 1: read 1
# Column 2: read 2 
# Column 3: output folder

# Submit job
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --machine-type n1-standard-2 \
  --preemptible \
  --logging gs://$OUTPUT_B/logging/ \
  --input R1="gs://$OUTPUT_B/A549-extract/split_fastq/A549-extract_L01.R1.001.fastq" \
  --input R2="gs://$OUTPUT_B/A549-extract/split_fastq/A549-extract_L01.R2.001.fastq" \
  --output FOLDER=gs://$OUTPUT_B/A549-extract/trimmed_fastq/* \
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
  --wait

########################## Align a pair of fastq shards ################################


# Submit job
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --machine-type n1-highcpu-16 \
  --disk-size 20 \
  --preemptible \
  --logging gs://$OUTPUT_B/logging/ \
  --input R1="gs://$OUTPUT_B/A549-extract/trimmed_fastq/A549-extract_L01.R1.001_val_1.fq" \
  --input R2="gs://$OUTPUT_B/A549-extract/trimmed_fastq/A549-extract_L01.R2.001_val_2.fq" \
  --input-recursive REF_GENOME="gs://$REF_DATA_B/grc37" \
  --output OUTPUT_DIR="gs://$OUTPUT_B/A549-extract/aligned_per_chard/*" \
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
  --wait


########################## Split chard's BAM by chromosome ################################

# Submit job
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --preemptible \
  --disk-size 50 \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --input BAM_FILES="gs://$OUTPUT_B/A549-extract/aligned_per_chard/A549-extract_L01.R1.001_val_1_bismark_bt2_pe.bam" \
  --output OUTPUT_DIR="gs://$OUTPUT_B/A549-extract/bam_per_chard_and_chr/*" \
  --script $HOME/GITHUB_REPOS/wgbs-asm/split_bam.sh \
  --wait


########################## Merge all BAMs by chromosome, clean them ################################

# Submit job
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --preemptible \
  --machine-type n1-highmem-8 \
  --disk-size 50 \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --env SAMPLE="A549-extract" \
  --env CHR="1" \
  --input BAM_FILES="gs://$OUTPUT_B/A549-extract/bam_per_chard_and_chr/*chr1.bam" \
  --output OUTPUT_DIR="gs://$OUTPUT_B/A549-extract/bam_per_chr/*" \
  --script $HOME/GITHUB_REPOS/wgbs-asm/merge_bam.sh \
  --wait


########################## Perform net methylation ################################

# Note: this step is not required to perform allele specific methylation

dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --preemptible \
  --machine-type n1-highcpu-16 \
  --disk-size 50 \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --input BAM="gs://$OUTPUT_B/A549-extract/bam_per_chr/A549-extract_chr1.bam" \
  --output OUTPUT_DIR="gs://$OUTPUT_B/A549-extract/net_methyl/*" \
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
  --wait


# Get some bedgraphs of net methylation and parameters through BigQuery


########################## SNP calling ################################

# Re-calibrate the BAM files.
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --preemptible \
  --machine-type n1-standard-16 \
  --disk-size 50 \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --env SAMPLE="A549-extract" \
  --env CHR="1" \
  --input BAM="gs://$OUTPUT_B/A549-extract/bam_per_chr/A549-extract_chr1_merged.bam" \
  --input REF_GENOME="gs://$REF_DATA_B/grc37"
  --input VCF="gs://$REF_DATA_B/dbSNP150_grc37_GATK/no_chr_dbSNP150_GRCh37.vcf"
  --output OUTPUT_DIR="gs://$OUTPUT_B/A549-extract/bam_per_chr/*" \
  --script bam_recalibration.sh \
  --wait


# Perform variant call.
