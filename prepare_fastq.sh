## This is an attempt at using the gcloud alpha genomics pipeline


# Create a bucket with the zipped fastq files from your dataset

# Create a directory called "logging"

pip install --upgrade google-api-python-client

https://cloud.google.com/genomics/docs/reference/gcloud-examples

https://codelabs.developers.google.com/codelabs/hpc-with-pipelines/index.html?index=..%2F..index#2


if [ "$FACILITY" == "encode" ] ; then
    while read SAMPLE ; do
        echo "Creating the file of lanes for the sample" $SAMPLE
        cat $HOME/encode_lanes.txt | grep -w $SAMPLE | \
        awk '!seen[$4]++' | awk '{ print $4 }' > $JOB_FOLDER/${SAMPLE}_lanes.List
    done < $JOB_FOLDER/sample.List
fi


SIZE_DRIVE="50"

gcloud config list
gcloud config set compute/region ""


gcloud alpha genomics pipelines run \
    --project hackensack-tyco \
    --command-line 'gunzip -c $ZIP > $UNZIP' \
    --preemptible \
    --cpus 2 \
    --boot-disk-size="50" \
    --inputs ZIP=gs://encode-wgbs/k562-test-pipeline/chunk1.fastq.gz \
    --outputs UNZIP=gs://encode-serverless/k562/chunk1.fastq \
    --logging gs://em-scratch/unzip.log


###########################################

gcloud alpha genomics pipelines run \
    --project hackensack-tyco \
    --docker-image "grc.io/hackensack-tyco/wgbs-asm" \
    --command-line 'gunzip -c $ZIP > $UNZIP' \
    --preemptible \
    --cpus 2 \
    --boot-disk-size="50" \
    --inputs ZIP=gs://encode-wgbs/k562-test-pipeline/chunk1.fastq.gz \
    --outputs UNZIP=gs://encode-serverless/k562/chunk1.fastq \
    --logging gs://em-scratch/unzip.log

gcloud alpha genomics pipelines run \
    --project hackensack-tyco \
    --docker-image "grc.io/hackensack-tyco/wgbs-asm" \
    --command-line 'samtools index ${BAM} ${BAI}' \
    --inputs BAM=gs://encode-serverless/A549/A549_split_aligned_L01.000.R1_val_1_bismark_bt2_pe_sorted.bam \
    --outputs BAI=gs://encode-serverless/A549/test.bam.bai



##############################################

gcloud alpha genomics pipelines run \
    --regions us-central1 \
    --logging gs://em-scratch/example.log \
    --command-line 'echo "${MESSAGE}"' \
    --inputs MESSAGE='Hello'



# Not very useful.
    --logging gs://encode-serverless/logging/pipeline_$(date +%Y_%m_%d___%H_%M_%S).log
    --logging gs://encode-serverless/logging/pipeline_$(date +%Y%m%d_%H%M%S).log \


gcloud alpha genomics pipelines run \
    --regions us-central1 \
    --logging gs://encode-serverless/example.log \
    --command-line 'echo "Hello, world!"'




############



## Not working

PYTHONPATH=/Users/emmanuel/GITHUB_REPOS/pipelines-api-examples python /Users/emmanuel/GITHUB_REPOS/pipelines-api-examples/compress/run_compress.py \
  --project hackensack-tyco \
  --zones "us-*" \
  --disk-size 200 \
  --operation "gunzip" \
  --input gs://encode-wgbs/k562/ENCFF336KJH.fastq.gz \
  --output s//encode-serverless/k562/ \
  --logging gs://encode-serverless/logging \
  --poll-interval 20



###################

# Using dsub

https://cloud.google.com/genomics/docs/tutorials/dsub



git clone https://github.com/googlegenomics/dsub.git
cd dsub

python setup.py install




#########################################################
############################################################################


#!/bin/bash
# Parameters to replace:
GOOGLE_CLOUD_PROJECT=GOOGLE_CLOUD_PROJECT
INPUT_PATTERN=gs://BUCKET/*.vcf
OUTPUT_TABLE=GOOGLE_CLOUD_PROJECT:BIGQUERY_DATASET.BIGQUERY_TABLE
TEMP_LOCATION=gs://BUCKET/temp

COMMAND="/opt/gcp_variant_transforms/bin/vcf_to_bq \
  --project ${GOOGLE_CLOUD_PROJECT} \
  --input_pattern ${INPUT_PATTERN} \
  --output_table ${OUTPUT_TABLE} \
  --temp_location ${TEMP_LOCATION} \
  --job_name vcf-to-bigquery \
  --runner DataflowRunner"

gcloud alpha genomics pipelines run \
  --project "${GOOGLE_CLOUD_PROJECT}" \
  --logging "${TEMP_LOCATION}/runner_logs_$(date +%Y%m%d_%H%M%S).log" \
  --zones us-west1-b \
  --service-account-scopes https://www.googleapis.com/auth/cloud-platform \
  --docker-image gcr.io/gcp-variant-transforms/gcp-variant-transforms \
  --command-line "${COMMAND}"



###########################################################################
###########################################################################


# spinning up tasks
echo "spinning up tasks"


rm -f operations
touch operations

gcloud alpha genomics pipelines run \
    --pipeline-file unzip-fastq.yaml \
    --logging gs://em-scratch/unzip.log \
    --inputs ZIP="gs://encode-wgbs/k562-test-pipeline/chunk1.fastq.gz", \
             UNZIP="gs://encode-serverless/k562/chunk1.fastq" 


echo -e "\ncomplete"

echo "waiting for all jobs to complete"

for op in `cat operations | cut -d / -f 2 | cut -d ']' -f 1`; do
      echo -n "."
      CMD="gcloud --format='value(done)' alpha genomics operations describe $op"
      while [[ $(eval ${CMD}) != "True" ]]; do echo -n "X"; sleep 5; done
done

echo -e "\nall jobs done"

echo "outputs:"
for x in {0..50}; do
      gsutil cat gs://${BUCKET_ID}/primes/logs/primes$x-stdout.log
done