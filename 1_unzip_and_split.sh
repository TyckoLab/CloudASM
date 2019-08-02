#!/bin/bash

# Define working directory and create it
WD="/data"
mkdir -p $WD

# Subdirectories
mkdir -p $WD/"single" | mkdir -p $WD/"split"

#
printf "\n Downloading the zipped fastq file... \n"
    time gsutil -o \
        'GSUtil:parallel_thread_count=1' \
        -o 'GSUtil:sliced_object_download_max_components=8' \
        cp $FASTQ $WD/"single"

#
printf "\n Downloading the zipped fastq file... \n"
gunzip $WD/"single"/*fastq.gz > $WD/"single"/$LANE.${READ}.fastq

#
printf "\n Splitting the fastq file in 12M chunks... \n"
split -l 12000000 --numeric-suffixes \
    --suffix-length=3 --additional-suffix=.${READ}.fastq\
    $WD/"single"/$LANE.${READ}.fastq \
    $WD/"split"/${LANE}.

# 
printf "\n Uploading the 12M chunks... \n"
gsutil -m \
    -o GSUtil:parallel_composite_upload_threshold=150M \
    cp $WD/"split"/* \
    gs://$BUCKET_UPLOAD/$SAMPLE/"split"/"fastq"
