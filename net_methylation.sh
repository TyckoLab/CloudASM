
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

# Path of where you downloaded the Github scripts
SCRIPTS="$HOME/GITHUB_REPOS/wgbs-asm/net_methylation"

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


# Create a dataset in big query named "wgbs_asm" (do this only once)
DATASET_ID="wgbs_asm"
bq --location=US mk -d \
    --description "Context files for net methylation" \
    $DATASET_ID


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

################################# Append context files (one OB and one OT per sample)

# Prepare TSV file
echo -e "--input CONTEXT_FILES\t--output CONTEXT_MERGED" > context.tsv

while read SAMPLE ; do
  echo -e "gs://$OUTPUT_B/${SAMPLE}/net_methyl/CpG_OT_${SAMPLE}_chr*.txt\tgs://$OUTPUT_B/${SAMPLE}/net_methyl/CpG_OT_${SAMPLE}_allchr.txt" >> context.tsv
  echo -e "gs://$OUTPUT_B/${SAMPLE}/net_methyl/CpG_OB_${SAMPLE}_chr*.txt\tgs://$OUTPUT_B/${SAMPLE}/net_methyl/CpG_OB_${SAMPLE}_allchr.txt" >> context.tsv
done < sample_id.txt

# Append all CpG_OB files
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --preemptible \
  --machine-type n1-standard-2 \
  --disk-size 50 \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --script ${SCRIPTS}/merge_context.sh \
  --tasks context.tsv \
  --wait


################################# Export context files in Big Query


# Prepare TSV file
echo -e "--env DATASET_ID\t--env SAMPLE\t--env STRAND\t--input CONTEXT" > context_to_bq.tsv

# Note: we replace the env variable with the same without any "dash"
while read SAMPLE ; do
  echo -e "${DATASET_ID}\t${SAMPLE}\tOT\tgs://$OUTPUT_B/${SAMPLE}/net_methyl/CpG_OT_${SAMPLE}_allchr.txt" >> context_to_bq.tsv
  echo -e "${DATASET_ID}\t${SAMPLE}\tOB\tgs://$OUTPUT_B/${SAMPLE}/net_methyl/CpG_OB_${SAMPLE}_allchr.txt" >> context_to_bq.tsv
done < sample_id.txt

# Export all (OB, OT) pairs per sample in Big Query
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --preemptible \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --command 'bq --location=US load \
               --replace \
               --source_format=CSV \
               --field_delimiter " " \
               ${DATASET_ID}.${SAMPLE}_CpG${STRAND} \
               ${CONTEXT} \
               read_id:STRING,meth_state:STRING,chr:STRING,pos:INTEGER,meth_call:STRING' \
  --tasks context_to_bq.tsv \
  --wait

## Note: the job fails if the field delimiter is the tab.

################################# Create bedgraph files of CpG coverage and CpG methylation % ########

# Prepare TSV file
echo -e "--env SAMPLE" > bedgraph.tsv

# Note: we replace the env variable with the same without any "dash"
while read SAMPLE ; do
  echo -e "${SAMPLE}" >> bedgraph.tsv
done < sample_id.txt

dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --preemptible \
  --zones $ZONE_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --env PROJECT_ID="${PROJECT_ID}" \
  --env DATASET_ID="${DATASET_ID}" \
  --env OUTPUT_B="$OUTPUT_B" \
  --script ${SCRIPTS}/bedgraph.sh \
  --tasks bedgraph.tsv \
  --wait


########### A few interesting notes

--gm12878
--590,000 CpGs in chr 22.
--528,800 with at least a coverage of 10x
--404,000 with at least a cov of 10x and a methylation of less than 100% and more than zero
-- 70% of CpG are left after these simple filters
