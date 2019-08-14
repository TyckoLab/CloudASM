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
       CAST(pos as INT64) - 250 as inf,
       CAST(pos as INT64) - 250 as sup,
       snp_id,
       ref,
       alt,
       CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) as ref_n,
       CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) as alt_n,
       CONCAT(chr,':', pos,'-',pos) as coord
    FROM ${DATASET_ID}.${SAMPLE}_chr${CHR}_vcf
    WHERE 
        snp_id LIKE '%rs%'
        AND CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) >= 10"

# Export variant list to the bucket
bq extract \
    --field_delimiter "\t" \
    --print_header=false \
    ${DATASET_ID}.${SAMPLE}_${CHR}_500bp \
    gs://$BUCKET/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}_500bp.bed
