#!/bin/bash

# Remove the rows where ref and alt are not single nucleotides
# Remove the rows where the snp_id is not in the form of "rs"

# The raw VCF file is a file used to remove CpGs that may overlap with SNPs identified in there
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
        ${DATASET_ID}.${SAMPLE}_vcf_raw_uploaded
    WHERE
        -- demand that the SNP is like rs-
        snp_id LIKE '%rs%'
        -- demand that we're deadling with a single nucleotide in REF and ALT columns
        AND BYTE_LENGTH(ref) = 1
        AND BYTE_LENGTH(alt) = 1
    "

# Delete the raw VCF table.
bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_raw_uploaded