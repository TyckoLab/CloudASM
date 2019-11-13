#!/bin/bash

# Create a table with a variant near a CpG in a 500bp window
# Remove the rows where ref and alt are not single nucleotides
# Remove the rows where the snp_id is not in the form of "rs"
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_vcf \
    --replace=true \
    "WITH
      -- We create a file with a 500bp window around the SNP and calculate the cov of the SNP
      variants AS (
        SELECT 
          chr,
          snp_id,
          ref,
          alt,
          SAFE_CAST(pos as INT64) as pos
        FROM 
          ${DATASET_ID}.${SAMPLE}_vcf_uploaded
        WHERE
          -- demand that the SNP is like rs-
          snp_id LIKE '%rs%'
          -- demand that we're deadling with a single nucleotide in REF and ALT columns
          AND BYTE_LENGTH(ref) = 1
          AND BYTE_LENGTH(alt) = 1
          -- below, we extract the coverage and demand that is it at least 10
        )
  -- Find the strand where the SNP can be identified in bisulfite-converted reads
  SELECT DISTINCT 
    snp_id,
    chr,
    pos,
    ref,
    alt,
    -- below, we indicate on which genome strand the SNP REF/ALT could be observed 
    -- with no ambiguity in a bisulfite converted sequence
    IF (ref = 'G' AND alt = 'A', TRUE,
            IF (ref = 'A' AND alt = 'G', TRUE,
              IF (ref = 'C' AND alt = 'T', FALSE,
                IF (ref = 'T' AND alt = 'C', FALSE, TRUE)))) AS CT_strand,
    IF (ref = 'G' AND alt = 'A', FALSE,
            IF (ref = 'A' AND alt = 'G', FALSE,
              IF (ref = 'C' AND alt = 'T', TRUE,
                IF (ref = 'T' AND alt = 'C', TRUE, TRUE)))) AS GA_strand 
  FROM
     variants
  "

# Delete VCF file that was uploaded
bq rm -f -t ${DATASET_ID}.${SAMPLE}_vcf_uploaded
