#!/bin/bash

# Build a table where each combination of read and variant where the variant is overlapping the read
# This will discard all paired reads not overlapping the variant but for which we know the genotyping

# Extract from the SAM file
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_chr${CHR}_${INF}_${SUP}_sam \
    --replace=true \
    " 
      SELECT 
        read_id,
        read_start,
        read_end,
        CT_strand AS seq_CT_strand,
        GA_strand AS seq_GA_strand,
        cigar,
        seq,
        score_before_recal
      FROM 
        ${DATASET_ID}.${SAMPLE}_recal_sam
      WHERE
        chr = '${CHR}'
        AND read_end <= ${SUP}
        AND read_start >= ${INF}
    "

# Extract from the VCF file
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_chr${CHR}_${INF}_${SUP}_vcf \
    --replace=true \
    "
      SELECT 
          snp_id,
          chr,
          pos,
          ref,
          alt,
          CT_strand,
          GA_strand
        FROM 
          ${DATASET_ID}.${SAMPLE}_vcf
        WHERE 
            chr = '${CHR}'
            AND pos <= ${SUP}
            AND pos >= ${INF}
    "

# Perform the JOIN
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_chr${CHR}_${INF}_${SUP}_tmp \
    --replace=true \
    "
      -- Find variants on DNA fragments where R1 or R2 has an overlap
      WITH sequences AS (SELECT * FROM ${DATASET_ID}.${SAMPLE}_chr${CHR}_${INF}_${SUP}_sam ),
           variants AS (SELECT * FROM ${DATASET_ID}.${SAMPLE}_chr${CHR}_${INF}_${SUP}_vcf )
        SELECT 
          snp_id,
          chr,
          pos,
          ref,
          alt,
          CT_strand,
          GA_strand,
          read_start,
          read_end,
          seq_CT_strand,
          seq_GA_strand,
          cigar,
          read_id,
          seq,
          score_before_recal
        FROM sequences
        INNER JOIN variants
        ON 
            read_start <= pos
            AND read_end >= pos
            AND (seq_CT_strand = CT_strand OR seq_GA_strand = GA_strand)
      "

# Delete intermediary tables
bq rm -f -t ${DATASET_ID}.${SAMPLE}_chr${CHR}_${INF}_${SUP}_sam
bq rm -f -t ${DATASET_ID}.${SAMPLE}_chr${CHR}_${INF}_${SUP}_vcf