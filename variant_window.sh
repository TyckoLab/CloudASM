#!/bin/bash

# Import CpG positions (only once)


# Import the VCF file in Big Query without the VCF header
bq --location=US load \
               --replace \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 117 \
               ${DATASET_ID}.${SAMPLE}_chr${CHR}_vcf \
               gs://$BUCKET/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}.vcf \
               chr:STRING,pos:STRING,snp_id:STRING,ref:STRING,alt:STRING,qual:FLOAT,filter:STRING,info:STRING,format:STRING,data:STRING

# Create a table with a min coverage of 10x per variant and a 500bp window
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_chr${CHR}_500bp_wCpG \
    --replace=true \
    "WITH
      win_2000 AS (
        SELECT 
          chr,
          snp_id,
          SAFE_CAST(pos AS INT64) - 250 AS inf_500, 
          SAFE_CAST(pos AS INT64) + 250 AS sup_500,
          SAFE_CAST(pos AS INT64) - 1000 AS inf_2000,
          SAFE_CAST(pos AS INT64) + 1000 AS sup_2000,
          ref,
          alt,
          SAFE_CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) as cov,
          SAFE_CAST(pos as INT64) as pos
        FROM 
          ${DATASET_ID}.${SAMPLE}_chr${CHR}_vcf    
        ),
      -- we convert the 2000bp window boundaries into string
      win_500 AS (
        SELECT 
          chr, 
          inf_500,
          sup_500,
          snp_id,
          ref,
          alt,
          cov,
          pos,
          SAFE_CAST(inf_2000 AS STRING) as inf_2000,
          SAFE_CAST(sup_2000 AS STRING) as sup_2000
        FROM win_2000
        ),
      -- create a coordinate of a 2000 bp window (to be used later)
      variants AS (
        SELECT 
          chr,
          inf_500,
          sup_500,
          snp_id,
          ref,
          alt,
          cov,
          pos,
          CONCAT(chr,':', inf_2000, '-', sup_2000) as coord_2000
        FROM win_500
        -- we ask that snp ID are rs* and that ref and alt nucleotides have only
        WHERE 
          snp_id LIKE '%rs%'
          AND cov >= 10
          AND BYTE_LENGTH(ref) = 1
          AND BYTE_LENGTH(alt) = 1
        ),
      cpg_pos AS (
        SELECT
          chr AS chr_cpg_pos,
          pos AS inf_cpg_pos,
          pos AS sup_cpg_pos
        FROM
          ${DATASET_ID}.${SAMPLE}_both_context
        GROUP BY
          chr, pos
        )
  -- We make sure that there is a CpG in the 500bp window of the SNP
  SELECT DISTINCT
    snp_id,
    chr,
    pos,
    coord_2000,
    ref,
    alt,
    cov,
    -- below, we indicate on which strand the SNP REF/ALT could be observed
    IF (ref = 'G' AND alt = 'A', 'CT',
            IF (ref = 'A' AND alt = 'G','CT',
              IF (ref = 'C' AND alt = 'T', 'GA',
                IF (ref = 'T' AND alt = 'C', 'GA', 'both')))) AS strand 
  FROM
     variants
  INNER JOIN
    cpg_pos ON
      cpg_pos.chr_cpg_pos = chr
      AND cpg_pos.inf_cpg_pos >= inf_500
      AND cpg_pos.sup_cpg_pos <= sup_500
  "


# Export variant list to the bucket
bq extract \
    --field_delimiter "\t" \
    --print_header=false \
    ${DATASET_ID}.${SAMPLE}_chr${CHR}_500bp_wCpG \
    gs://$BUCKET/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}_variants.txt
