
########################## Copy and paste within the bash ################################

# GCP global variables
PROJECT_ID="hackensack-tyco"
REGION_ID="us-central1"
ZONE_ID="us-central1-b"
DOCKER_IMAGE="gcr.io/hackensack-tyco/wgbs-asm"

# Big Query variables
DATASET_ID="wgbs_asm" 

# Cloud storage variables
OUTPUT_B="em-encode-deux" # will be created by the script
REF_DATA_B="wgbs-ref-files" # See documentation for what it needs to contain

# Path of where you downloaded the Github scripts
SCRIPTS="$HOME/GITHUB_REPOS/wgbs-asm/"

# Create a bucket with the analysis
gsutil mb -c regional -l $REGION_ID gs://$OUTPUT_B 

# Create a local directory to download files
WD="$HOME/wgbs" && mkdir -p $WD && cd $WD

# Download the meta information about the samples and files to be analyzed.
gsutil cp gs://$INPUT_B/samples.tsv $WD
dos2unix samples.tsv 

# List of samples
awk -F "\t" \
    '{if (NR!=1) \
    print $1}' samples.tsv | uniq > sample_id.txt

echo "There are" $(cat sample_id.txt | wc -l) "samples to be analyzed"

########################## Unzip, rename, and split fastq files ################################

# Create an TSV file with parameters for the job
awk -v INPUT_B="${INPUT_B}" \
    -v OUTPUT_B="${OUTPUT_B}" \
    'BEGIN { FS=OFS="\t" } 
    {if (NR!=1) 
        print "gs://"INPUT_B"/"$1"/"$2, 
              $5".fastq", 
              "gs://"OUTPUT_B"/"$1"/split_fastq/*.fastq" 
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
             split -l 1200000 \
                --numeric-suffixes --suffix-length=4 \
                --additional-suffix=.fastq \
                $(dirname "${ZIPPED}")/${FASTQ} \
                $(dirname "${OUTPUT_FILES}")/${FASTQ%fastq}' \
  --tasks decompress.tsv \
  --wait

########################## Trim a pair of fastq shards ################################

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

# Submit job. Make sure you have enough resources
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
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

# Submit job
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --machine-type n1-standard-16 \
  --preemptible \
  --disk-size 40 \
  --logging gs://$OUTPUT_B/logging/ \
  --input-recursive REF_GENOME="gs://$REF_DATA_B/grc37" \
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
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --script ${SCRIPTS}/split_bam.sh \
  --tasks split_bam.tsv \
  --wait


########################## Merge all BAMs by chromosome, clean them ################################

# Prepare TSV file
echo -e "--env SAMPLE\t--env CHR\t--input BAM_FILES\t--output OUTPUT_DIR" > merge_bam.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
  echo -e "${SAMPLE}\t${CHR}\tgs://$OUTPUT_B/$SAMPLE/bam_per_chard_and_chr/*chr${CHR}.bam\tgs://$OUTPUT_B/$SAMPLE/bam_per_chr/*" >> merge_bam.tsv
  done
done < sample_id.txt

# Print a message in the terminal
echo "There are" $(cat merge_bam.tsv | wc -l) "to be launched"

# Submit job
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --machine-type n1-highmem-8 \
  --preemptible \
  --disk-size 30 \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --script ${SCRIPTS}/merge_bam.sh \
  --tasks merge_bam.tsv \
  --wait

################################# Net methylation call

# Note: this step requires one BAM file per chromosome (in "bam_per_chr") 

# Prepare TSV file
echo -e "--input BAM\t--output OUTPUT_DIR" > methyl.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
  echo -e "gs://$OUTPUT_B/$SAMPLE/bam_per_chr/${SAMPLE}_chr${CHR}.bam\tgs://$OUTPUT_B/$SAMPLE/net_methyl/*" >> methyl.tsv
  done
done < sample_id.txt

# Print a message in the terminal
echo "There are" $(tail -n +2 methyl.tsv | wc -l) "to be launched"

dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --preemptible \
  --machine-type n1-highmem-8 \
  --disk-size 50 \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
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
  --preemptible \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --command 'bq --location=US load \
               --replace=false \
               --source_format=CSV \
               --skip_leading_rows 1 \
               --field_delimiter " " \
               ${DATASET_ID}.${SAMPLE}_CpG${STRAND} \
               ${CONTEXT} \
               read_id:STRING,meth_state:STRING,chr:STRING,pos:INTEGER,meth_call:STRING' \
  --tasks context_to_bq.tsv \
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
  --image $DOCKER_IMAGE \
  --preemptible \
  --logging gs://$OUTPUT_B/logging/ \
  --env PROJECT_ID="${PROJECT_ID}" \
  --env DATASET_ID="${DATASET_ID}" \
  --env OUTPUT_B="$OUTPUT_B" \
  --script ${SCRIPTS}/append_context.sh \
  --tasks append_context.tsv \
  --wait


########################## Re-calibrate BAM  ################################

# This step is required by the variant call Bis-SNP

# Prepare TSV file
echo -e "--env SAMPLE\t--env CHR\t--input BAM\t--output OUTPUT_DIR" > bam_recalibration.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
  echo -e "$SAMPLE\t$CHR\tgs://$OUTPUT_B/$SAMPLE/bam_per_chr/${SAMPLE}_chr${CHR}.bam\tgs://$OUTPUT_B/$SAMPLE/recal_bam_per_chr/*" >> bam_recalibration.tsv
  done
done < sample_id.txt

# Print a message in the terminal
echo "There are" $(tail -n +2 bam_recalibration.tsv | wc -l) "to be launched"

# Re-calibrate the BAM files.
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --machine-type n1-standard-16 \
  --disk-size 200 \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --input REF_GENOME="gs://$REF_DATA_B/grc37/*" \
  --input VCF="gs://$REF_DATA_B/dbSNP150_grc37_GATK/no_chr_dbSNP150_GRCh37.vcf" \
  --script ${SCRIPTS}/bam_recalibration.sh \
  --tasks bam_recalibration.tsv \
  --wait


########################## Variant call  ################################

# Prepare TSV file
echo -e "--env SAMPLE\t--env CHR\t--input BAM_BAI\t--output OUTPUT_DIR" > variant_call.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
  echo -e "$SAMPLE\t$CHR\tgs://$OUTPUT_B/${SAMPLE}/recal_bam_per_chr/${SAMPLE}_chr${CHR}_recal.ba*\tgs://$OUTPUT_B/${SAMPLE}/variants_per_chr/*" >> variant_call.tsv
  done
done < sample_id.txt


dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --machine-type n1-standard-16 \
  --disk-size 500 \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --input REF_GENOME="gs://$REF_DATA_B/grc37/*" \
  --input VCF="gs://$REF_DATA_B/dbSNP150_grc37_GATK/no_chr_dbSNP150_GRCh37.vcf" \
  --script ${SCRIPTS}/variant_call.sh \
  --tasks variant_call.tsv \
  --wait


########################## TEMP TO BE DELETED ################################

#### TEMP TO BE DELETED

# dsub \
#   --provider google-v2 \
#   --project $PROJECT_ID \
#   --preemptible \
#   --machine-type n1-standard-4 \
#   --disk-size 500 \
#   --zones $ZONE_ID \
#   --image $DOCKER_IMAGE \
#   --logging gs://$OUTPUT_B/logging/ \
#   --input BAM="gs://$OUTPUT_B/${SAMPLE}/recal_bam_per_chr/${SAMPLE}_chr${CHR}_recal.bam" \
#   --input BAI="gs://$OUTPUT_B/${SAMPLE}/recal_bam_per_chr/${SAMPLE}_chr${CHR}_recal.bai" \
#   --output SAM="gs://$OUTPUT_B/${SAMPLE}/recal_bam_per_chr/${SAMPLE}_chr${CHR}_recal.sam" \
#   --command 'samtools view -o \
#                 ${SAM} \
#                 ${BAM}' \
#   --wait
##############


########################## Export recal bam to Big Query ################################

# The SAM file was created by the variant_call script

# Prepare TSV file
echo -e "--env SAMPLE\t--env SAM" > sam_to_bq.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
    echo -e "$SAMPLE\tgs://$OUTPUT_B/${SAMPLE}/recal_bam_per_chr/${SAMPLE}_chr${CHR}_recal.sam" >> sam_to_bq.tsv
  done
  
  # Delete existing SAM on big query
  bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_recal_sam_raw

done < sample_id.txt

# We append all chromosomes in the same file.
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --command 'bq --location=US load \
               --replace=false \
               --source_format=CSV \
               --field_delimiter "\t" \
               ${DATASET_ID}.${SAMPLE}_recal_sam_raw \
               ${SAM} \
               read_id:STRING,flag:INTEGER,chr:STRING,read_start:INTEGER,mapq:INTEGER,cigar:STRING,rnext:STRING,mate_read_start:INTEGER,insert_length:INTEGER,seq:STRING,score:STRING,bismark:STRING,picard_flag:STRING,read_g:STRING,genome_strand:STRING,NM_tag:STRING,meth:STRING,score_before_recal:STRING,read_strand:STRING' \
  --tasks sam_to_bq.tsv \
  --wait

########################## Clean each SAM (one per sample) %=##################

# Prepare TSV file
echo -e "--env SAMPLE" > clean_sam.tsv

while read SAMPLE ; do
  echo -e "$SAMPLE" >> clean_sam.tsv
done < sample_id.txt


dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --env PROJECT_ID="${PROJECT_ID}" \
  --script ${SCRIPTS}/clean_sam.sh \
  --tasks clean_sam.tsv \
  --wait


# Delete the SAM files from the bucket (they take a lot of space)
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --preemptible \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --command 'gsutil rm ${SAM}' \
  --tasks sam_to_bq.tsv \
  --wait


########################## Export VCF to Big Query (one per sample) ################################


# Prepare TSV file
echo -e "--env SAMPLE\t--env VCF" > vcf_to_bq.tsv

while read SAMPLE ; do
  for CHR in `seq 1 22` X Y ; do 
    echo -e "$SAMPLE\tgs://$BUCKET/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}.vcf" >> vcf_to_bq.tsv
  done
  
  # Delete existing SAM on big query
  bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_raw

done < sample_id.txt

# We append all chromosomes files in the same file.
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --command 'bq --location=US load \
               --replace=false \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 117 \
               ${DATASET_ID}.${SAMPLE}_vcf_raw \
               ${VCF} \
               chr:STRING,pos:STRING,snp_id:STRING,ref:STRING,alt:STRING,qual:FLOAT,filter:STRING,info:STRING,format:STRING,data:STRING' \
  --tasks vcf_to_bq.tsv \
  --wait

########################## Clean each VCF (one per sample) %=##################

# Filter out the SNPs that are not within 500bp of a CpG that is at least 10x covered.
# This removes about 5% of SNPs. Takes ~30min

# REMOVE THE INNER JOIN WITH CPG SITES AS IT TAKES A LONG TIME TO REMOVE JUST 5 PC

# Prepare TSV file
echo -e "--env SAMPLE" > clean_vcf.tsv

while read SAMPLE ; do
  echo -e "$SAMPLE" >> clean_vcf.tsv
done < sample_id.txt


dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --env PROJECT_ID="${PROJECT_ID}" \
  --script ${SCRIPTS}/clean_vcf.sh \
  --tasks clean_vcf.tsv \
  --wait


########################## Find the read IDs that overlap the snp ##################

# Prepare TSV file
echo -e "--env SAMPLE\t--env CHR" > reads_overlap_snp.tsv

while read SAMPLE ; do
  for CHR in `seq 21 22`  ; do 
    echo -e "${SAMPLE}\t${CHR}" >> reads_overlap_snp.tsv
  done
  
  # Delete existing SAM on big query
  bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_reads  

done < sample_id.txt

# Create one file per chromosome
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --env PROJECT_ID="${PROJECT_ID}" \
  --script ${SCRIPTS}/reads_overlap_snp.sh \
  --tasks reads_overlap_snp.tsv \
  --wait

# Merge all chromosome files into a single file per sample 
# and delete the chromosome file
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --env DATASET_ID="${DATASET_ID}" \
  --env PROJECT_ID="${PROJECT_ID}" \
  --command 'bq cp \
               --append_table \
               ${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_${CHR} \
               ${DATASET_ID}.${SAMPLE}_vcf_reads \
             && bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_${CHR}' \
  --tasks reads_overlap_snp.tsv \
  --wait



########################## Tag each read with REF or ALT ##################

# We consider the cases where there are 1, 3, and 5 numbers in the CIGAR string
# We leave out the 0.00093% where the CIGAR string has 7 numbers or more
# Note: 0.2% of SNPs are left out.


96% are a perfect match
they all start with M
2.5% of (snp,read_id) need to deal with the CIGAR

########################## Split the variants in 200 shards ################################

# We use the 500bp window created above to qualify for "near"

# Prepare TSV file
# echo -e "--env SAMPLE\t--env CHR\t--input VARIANTS_CHR\t--output OUTPUT_DIR" > variant_list.tsv

# while read SAMPLE ; do
#   for CHR in `seq 1 22` X Y ; do 
#   echo -e "${SAMPLE}\t${CHR}\tgs://${OUTPUT_B}/${SAMPLE}/variants_per_chr/${SAMPLE}_chr${CHR}_variants.txt\tgs://$OUTPUT_B/${SAMPLE}/variant_shards/*" >> variant_list.tsv
#   done
# done < sample_id.txt

# dsub \
#   --provider google-v2 \
#   --project $PROJECT_ID \
#   --zones $ZONE_ID \
#   --image $DOCKER_IMAGE \
#   --logging gs://$OUTPUT_B/logging/ \
#   --command 'split -l 200 \
#                 --numeric-suffixes --suffix-length=6 \
#                 --additional-suffix=.txt \
#                 ${VARIANTS_CHR}\
#                 $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_variants_' \
#   --tasks variant_list.tsv \
#   --wait

########################## Create a pair of (REF, ALT) of each BAM for each combination [SAM, snp shard] ################################

# Need SNP ID, chr, position, window of 2000 centered on SNP.
#100,000 SNP x 54 in chr 1

# SAMPLE="gm12878"
# CHR="22"

# dsub \
#   --provider google-v2 \
#   --project $PROJECT_ID \
#   --zones $ZONE_ID \
#   --machine-type n1-standard-4 \
#   --disk-size 20 \
#   --image $DOCKER_IMAGE \
#   --logging gs://$OUTPUT_B/logging/ \
#   --input BAM_BAI="gs://$OUTPUT_B/${SAMPLE}/recal_bam_per_chr/${SAMPLE}_chr${CHR}_recal.ba*" \
#   --input SNP_LIST="gs://$OUTPUT_B/$SAMPLE/variant_shards/split-1-of-50_${SAMPLE}_chr${CHR}.txt" \
#   --input VCF="gs://$OUTPUT_B/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}.vcf" \
#   --output OUTPUT_FOLDER="gs://$OUTPUT_B/$SAMPLE/genotype/*" \
#   --script ${SCRIPTS}/genotype.sh \
#   --tasks genotype.tsv \
#   --wait

