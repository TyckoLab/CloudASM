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

# prepare a batch file for decompressing the fastq files
awk '{if (NR!=1) {print}}' samples.tsv > samples_noheader.tsv

# Create an input file with parameters ()
awk -v INPUT_B="${INPUT_B}" \
    -v OUTPUT_B="${OUTPUT_B}" \
    'BEGIN { FS="\t"; OFS="\t" } \
    {print "gs://"INPUT_B"/"$1"/"$2, \
            10 * $6, \
           "gs://"OUTPUT_B"/"$1"/"$5".fastq"}' \
    samples_noheader.tsv > decompress.tsv 

sed -i '1i --input ZIPPED\t--env DISK_SIZE\t--output UNZIPPED' decompress.tsv


dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --logging gs://$OUTPUT_B/logging/ \
  --disk-size 10 \
  --preemptible \
  --image $DOCKER_IMAGE \
  --command 'gunzip ${ZIPPED} && mv ${ZIPPED%.gz} ${UNZIPPED}' \
  --tasks decompress.tsv \
  --wait



