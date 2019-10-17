#!/bin/bash

# Delete any previous database on Big Query
bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_snps_for_cpg_raw 

if [ "${SNPS_FOR_CPG}" = "common_snp" ]; then
    if [ "${GENOME}" = "hg19" ]; then
        COMMON_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151Common.txt.gz"
    else 
        echo "Downloading the common SNPs for GRCh38"
        COMMON_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz"
    fi

    # Download the common snps
    wget $COMMON_URL

    # Unzip the file
    gunzip $(basename "${COMMON_URL}")

    # Copy to bucket
    gsutil cp $(basename "${COMMON_URL%.gz}") gs://$OUTPUT_B/$SAMPLE/"snps_for_cpg"/$(basename "${COMMON_URL%.gz}")

    # Upload file to Big Query
    bq --location=US load \
        --replace=true \
        --source_format=CSV \
        --field_delimiter "\t" \
        ${DATASET_ID}.${SAMPLE}_snps_for_cpg_raw \
        gs://$OUTPUT_B/$SAMPLE/"snps_for_cpg"/$(basename "${COMMON_URL%.gz}") \
        bin:INTEGER,chr:STRING,pos:INTEGER,pos_sup:INTEGER,snp_id:STRING,score:INTEGER,strand:STRING,refNCBI:STRING,refUCSC:STRING,observed:STRING,molType:STRING,class:STRING,valid:STRING,avHet:FLOAT,avHetSE:FLOAT,func:STRING,locType:STRING,weight:INTEGER,exceptions:STRING,submitterCount:INTEGER,submitters:STRING,alleleFreqCount:INTEGER,alleles:STRING,alleleNs:STRING,alleleFreqs:STRING,bitfields:STRING


    # Clean the database of common snps
    bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_snps_for_cpg \
    --replace=true \
    "WITH
    PRE_CLEANING AS (
        SELECT
            pos,
            snp_id,
            alleleFreqs,
            REGEXP_REPLACE(chr, 'chr','') AS chr_clean,
            REGEXP_EXTRACT_ALL(alleleFreqs, r'[0-9]*\.?[0-9]+') AS allele
        FROM
            ${DATASET_ID}.${SAMPLE}_snps_for_cpg_raw
        WHERE
            snp_id LIKE '%rs%'
        ),
    CLEAN_DATA AS (
        SELECT
            REGEXP_REPLACE(chr_clean, r'_([a-zA-Z0-9\s]+)', '') AS chr,
            pos,
            snp_id,
            (SELECT MAX(x) AS max_freq
             FROM UNNEST(allele) AS x) AS max_freq
        FROM
            PRE_CLEANING
        ),
    FORMAT_DATA AS (
        SELECT
            chr AS chr_snp,
            pos AS pos_snp,
            snp_id,
            SAFE_CAST(max_freq AS FLOAT64) AS max
        FROM
            CLEAN_DATA 
        )
    SELECT *
    FROM FORMAT_DATA
    WHERE max < 1 - ${SNP_FREQ}
    "

else 
  echo "Using one of the Bis-SNP's VCF databases to remove CpGs"

  for CHR in `seq 1 22` X Y ; do 
    echo "Processing chr " $CHR
    sleep 1s
    bq --location=US load \
        --replace=false \
        --source_format=CSV \
        --field_delimiter "\t" \
        --skip_leading_rows 116 \
        ${DATASET_ID}.${SAMPLE}_snps_for_cpg_raw \
        gs://$OUTPUT_B/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}_${SNPS_FOR_CPG} \
        chr:STRING,pos:STRING,snp_id:STRING,ref:STRING,alt:STRING,qual:FLOAT,filter:STRING,info:STRING,format:STRING,data:STRING
  done

  # Clean the VCF of Bis-SNP
  bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_snps_for_cpg \
    --replace=true \
    "
    SELECT 
        chr AS chr_snp,
        SAFE_CAST(pos as INT64) AS pos_snp,
        snp_id,
        ref,
        alt
    FROM 
        ${DATASET_ID}.${SAMPLE}_snps_for_cpg_raw
    WHERE
        -- demand that the SNP is like rs-
        snp_id LIKE '%rs%'
        -- demand that we're deadling with a single nucleotide in REF and ALT columns
        AND BYTE_LENGTH(ref) = 1
        AND BYTE_LENGTH(alt) = 1
    "
fi

bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_snps_for_cpg_raw 
