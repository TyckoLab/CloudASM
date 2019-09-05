#!/bin/bash

# Create a table with a variant near a CpG in a 500bp window
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_chr${CHR}_tmp \
    --replace=true \
    "WITH
      -- We create a file with a 500bp window around the SNP and calculate the cov of the SNP
      variants AS (
        SELECT 
          chr,
          snp_id,
          SAFE_CAST(pos AS INT64) - 250 AS inf_500, 
          SAFE_CAST(pos AS INT64) + 250 AS sup_500,
          ref,
          alt,
          SAFE_CAST(pos as INT64) as pos
        FROM 
          ${DATASET_ID}.${SAMPLE}_vcf_raw
        WHERE
          -- demand that the SNP is like rs-
          snp_id LIKE '%rs%'
          -- demand that we're deadling with a single nucleotide in REF and ALT columns
          AND BYTE_LENGTH(ref) = 1
          AND BYTE_LENGTH(alt) = 1
          -- below, we extract the coverage and demand that is it at least 10
          AND chr = '${CHR}'
        ),
      cpg_pos AS (
        SELECT
          chr AS chr_cpg_pos,
          pos AS inf_cpg_pos,
          pos AS sup_cpg_pos
        FROM
          ${DATASET_ID}.${SAMPLE}_context
        WHERE
          chr = '${CHR}'
        GROUP BY
          chr, pos
        )
  -- We make sure that there is at least a variant that is 10x covered within 500bp of a CpG
  -- since there are several CpG that may qualify we keep the distinct results
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
  -- Making sure that there are CpG within 500bp of the SNP
  INNER JOIN
    cpg_pos ON
      cpg_pos.chr_cpg_pos = chr
      AND cpg_pos.inf_cpg_pos >= inf_500
      AND cpg_pos.sup_cpg_pos <= sup_500
  "

