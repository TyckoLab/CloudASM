#!/bin/bash

# Create a table with a min coverage of 10x per variant and a 500bp window
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf \
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
          SAFE_CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) as cov,
          SAFE_CAST(pos as INT64) as pos
        FROM 
          ${DATASET_ID}.${SAMPLE}_vcf_raw
        WHERE
          snp_id LIKE '%rs%'
          AND BYTE_LENGTH(ref) = 1
          AND BYTE_LENGTH(alt) = 1
          AND SAFE_CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) >= 10
        ),
      cpg_pos AS (
        SELECT
          chr AS chr_cpg_pos,
          pos AS inf_cpg_pos,
          pos AS sup_cpg_pos
        FROM
          ${DATASET_ID}.${SAMPLE}_context
        GROUP BY
          chr, pos
        )
  -- We make sure that there is at least a CpG that is 10x covered within 500bp of the SNP
  -- since there are several CpG that may qualify we keep the distinct results
  SELECT 
    DISTINCT snp_id,
    chr,
    pos,
    ref,
    alt,
    cov,
    -- below, we indicate on which genome strand the SNP REF/ALT could be observed 
    -- with no ambiguity in a bisulfite converted sequence
    IF (ref = 'G' AND alt = 'A', 'CT',
            IF (ref = 'A' AND alt = 'G','CT',
              IF (ref = 'C' AND alt = 'T', 'GA',
                IF (ref = 'T' AND alt = 'C', 'GA', 'both')))) AS genome_strand 
  FROM
     variants
  INNER JOIN
    cpg_pos ON
      cpg_pos.chr_cpg_pos = chr
      AND cpg_pos.inf_cpg_pos >= inf_500
      AND cpg_pos.sup_cpg_pos <= sup_500
  "

