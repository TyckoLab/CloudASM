#!/bin/bash

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
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_${CHR}_500bp \
    --replace=true \
    "SELECT 
       chr, 
       CAST(pos as INT64) - 250 as inf_500,
       CAST(pos as INT64) + 250 as sup_500,
       snp_id,
       ref,
       alt,
       CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) as ref_n,
       CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) as alt_n,
       CAST(pos as INT64) as pos,
       CONCAT(chr,':', CAST(pos as INT64) - 1000,'-',CAST(pos as INT64) + 1000)) as coord_2000
    FROM ${DATASET_ID}.${SAMPLE}_chr${CHR}_vcf
    WHERE 
        snp_id LIKE '%rs%'
        AND CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) >= 10"

### 
WITH win AS 
  (SELECT
    chr AS chr_win,
    snp_id AS snpid_win,
    CAST(pos as INT64) - 1000 as inf_2000,
    CAST(pos as INT64) + 1000 as sup_2000
    FROM `hackensack-tyco.wgbs_asm.gm12878_chr22_vcf`)
SELECT 
       chr, 
       CAST(pos as INT64) - 250 as inf_500,
       CAST(pos as INT64) + 250 as sup_500,
       snp_id,
       ref,
       alt,
       CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) as ref_n,
       CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) as alt_n,
       CAST(pos as INT64) as pos,
       CONCAT(chr,':', CAST((SELECT inf_2000 FROM win) AS STRING),'-', CAST((SELECT sup_2000 FROM win) AS STRING)) as coord_2000
       FROM `hackensack-tyco.wgbs_asm.gm12878_chr22_vcf`
UNION ALL 
WHERE 
  (SELECT snpid_win FROM win) = snp_id
  AND (SELECT chr_win FROM win) = chr
LIMIT 10

# Export variant list to the bucket
bq extract \
    --field_delimiter "\t" \
    --print_header=false \
    ${DATASET_ID}.${SAMPLE}_${CHR}_500bp \
    gs://$BUCKET/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}_500bp.bed
