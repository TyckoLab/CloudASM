# Pre-requisities


##### Data
# Have a bucket already created with one folder per sample. 
# Each folder should include the zipped fastq files
# In the bucket, there should be a tsv file "samples.tsv" with each file name, its lane number, its read number (R1 or R2), and the size of the zipped fastq in GB

##### Authentification
# Be authentificated by gcloud (type "gcloud config list" to confirm)

########################## Install a vritual env to run Python ################################

# Install Docker 
https://docs.docker.com/install/

# Install dsub

git clone https://github.com/googlegenomics/dsub.git
cd dsub

python setup.py install

# Install virtualenv, a tool to create isolated Python environments. 
# Install https://virtualenv.pypa.io/en/stable/installation/

# Go into any folder and type:
virtualenv --python=python2.7 dsub_libs

# Launch the virtual environment
source dsub_libs/bin/activate

########################## Copy and paste within the bash ################################

# GCP variables
PROJECT_ID="hackensack-tyco"
REGION_ID="us-central1"
ZONE_ID="us-central1-b"
DOCKER_IMAGE="gcr.io/hackensack-tyco/wgbs-asm"

# Names of the buckets (must be unique)
INPUT_B="encode-wgbs" # should already exist with the zipped fastq files inside
OUTPUT_B="em-encode" # will be created by the script

# Create a bucket with the analysis
gsutil mb -c coldline -l $REGION_ID gs://$OUTPUT_B 

# Create a local directory to download files
WD="$HOME/wgbs" && mkdir -p $WD && cd $WD

# Download the meta information about the samples and files to be analyzed.
gsutil cp gs://$INPUT_B/samples.tsv $WD

# Check the number of samples to be analyzed.
echo "There are" $(($(awk -F '\t' '{print $1}' samples.tsv | uniq -c | wc -l) - 1)) "samples to be analyzed"

# Find the largest fastq file and size disks accordingly.

########################## Simple test for dsub ################################

# This should not lead to any error
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --regions $REGION_ID \
  --logging gs://$OUTPUT_B/logging/ \
  --output OUT=gs://$OUTPUT_B/output/out.txt \
  --command 'echo "Hello World" > "${OUT}"' \
  --wait

# Delete the output
gsutil rm -r gs://$OUTPUT_B/output

########################## Unzip and rename fastq files ################################

# Create an input file with parameters ()
# awk -v INPUT_B="${INPUT_B}" \
#     -v OUTPUT_B="${OUTPUT_B}" \
#     'BEGIN { FS="\t"; OFS="\t" } \
#     {if (NR!=1) \
#      print $1, \
#            "gs://"INPUT_B"/"$1"/"$2, \
#            "gs://"OUTPUT_B"/"$1"/"$5".fastq"}' \
#     samples.tsv > decompress.tsv 

# # Add headers to the file
# sed -i '1i --env SAMPLE\t--input ZIPPED\t--output UNZIPPED' decompress.tsv

# Unzip the fastq files (WORKING)
# dsub \
#   --provider google-v2 \
#   --project $PROJECT_ID \
#   --zones $ZONE_ID \
#   --logging gs://$OUTPUT_B/logging/ \
#   --disk-size 10 \
#   --preemptible \
#   --image $DOCKER_IMAGE \
#   --command 'gunzip ${ZIPPED} && mv ${ZIPPED%.gz} ${UNZIPPED}' \
#   --tasks decompress.tsv \
#   --wait

awk -v INPUT_B="${INPUT_B}" \
    -v OUTPUT_B="${OUTPUT_B}" \
    'BEGIN { FS="\t"; OFS="\t" } \
    {if (NR!=1) \
     print "gs://"INPUT_B"/"$1"/"$2, \
           $5".fastq", \
           "gs://"OUTPUT_B"/"$1"/split/*.fastq"}' \
    samples.tsv > decompress.tsv 

# Add headers to the file
sed -i '1i --input ZIPPED\t--env FASTQ\t--output OUTPUT_FILES' decompress.tsv


# Unzip and split the fastq files in 12M-row files.
# NOT WORKING
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
                --numeric-suffixes --suffix-length=3 \
                --additional-suffix=.fastq \
                $(dirname "${ZIPPED}")/${FASTQ} \
                $(dirname "${OUTPUT_FILES}")/${FASTQ%fastq}' \
  --tasks decompress.tsv \
  --wait


dsub \
  --provider local \
  --logging $WD/example/logging/ \
  --image $DOCKER_IMAGE \
  --input ZIPPED=$WD/example/sample1_chunk1.fastq.gz \
  --output OUTPUT_FILES=$WD/example/split/*.fastq \
  --env FASTQ="sample1_L01.R1.fastq" \
  --command 'gunzip ${ZIPPED} && \
              mv ${ZIPPED%.gz} $(dirname "${ZIPPED}")/$FASTQ && \
              split -l 120000 \
                  --numeric-suffixes --suffix-length=3 \
                  --additional-suffix=.fastq \
                  $(dirname "${ZIPPED}")/${FASTQ} \
                  $(dirname "${OUTPUT_FILES}")/${FASTQ%fastq}' \
  --wait

########################## Split the fastq files ################################

# Fetch the names of the unzipped fastq files
awk -v OUTPUT_B="${OUTPUT_B}" \
    'BEGIN { FS="\t"; OFS="\t" } \
    {if (NR!=1) \
    print $3, \
          "gs://"OUTPUT_B"/"$1"/split/*.fastq"}' \
    decompress.tsv > split.tsv

# Add headers to the file
sed -i '1i --input FASTQ\t--output OUTPUT_FILES' split.tsv

# Split the fastq files in 12M-row files.
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --logging gs://$OUTPUT_B/logging/ \
  --disk-size 10 \
  --preemptible \
  --image $DOCKER_IMAGE \
  --command 'split -l 120000 \
             --numeric-suffixes --suffix-length=3 \
             --additional-suffix=.fastq \
             ${FASTQ} $(dirname "${OUTPUT_FILES}")/$(basename "${FASTQ%fastq}")' \
  --tasks split.tsv \
  --wait

# dsub \
#   --provider local \
#   --logging $WD/logging/ \
#   --image $DOCKER_IMAGE \
#   --input FASTQ=$WD/sample1_L01.R1.fastq \
#   --output OUTPUT_FILES=$WD/split/*.fastq \
#   --command 'split -l 120000 \
#              --numeric-suffixes --suffix-length=3 \
#              --additional-suffix=.fastq \
#              ${FASTQ} $(dirname "${OUTPUT_FILES}")/$(basename "${FASTQ%fastq}")' \
#   --wait


split -l 12000000 --numeric-suffixes \
    --suffix-length=3 --additional-suffix=.R2.fastq\
    $WD/"unzipped_fastq"/${LANE}.R2.fastq \
    $WD/"split_fastq"/${LANE}.
