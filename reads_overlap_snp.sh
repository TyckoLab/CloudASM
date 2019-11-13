#!/bin/bash

# First step: find which unique read_id (R1 and R2) each variant is overlapping with
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_chr${CHR}_first \
    --replace=true \
    " WITH
        sequences AS (
          SELECT 
            read_id,
            read_start,
            read_end,
            CT_strand AS seq_CT_strand,
            GA_strand AS seq_GA_strand
          FROM 
            ${DATASET_ID}.${SAMPLE}_recal_sam
          WHERE
            chr = '${CHR}'
        ),
       variants AS (
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
      ),
      -- Find the lower and upper bounds of a pair of read id (R1 and R2)
      unique_read_id AS (
        SELECT 
          read_id AS read_id_unique,
          min(read_start) AS seq_start,
          max(read_end) AS seq_end,
          ANY_VALUE(seq_CT_strand) AS seq_CT_strand,
          ANY_VALUE(seq_GA_strand) AS seq_GA_strand
        FROM sequences
        GROUP BY read_id
      )
      -- Find variants on DNA fragments where R1 or R2 has an overlap
        SELECT 
          snp_id, 
          pos, 
          ref, 
          alt, 
          read_id_unique, 
          CT_strand, 
          GA_strand 
        FROM unique_read_id
        INNER JOIN variants
        ON 
            seq_start <= pos
            AND seq_end >= pos
            AND (seq_CT_strand = CT_strand OR seq_GA_strand = GA_strand)
      "

# Second step: associate these variants found above with the read id that do not overlap with them (but overlap with their respective R1 or R2)
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_chr${CHR} \
    --replace=true \
      "
      WITH
        sequences AS (
          SELECT 
            read_id,
            read_start,
            read_end,
            chr,
            CT_strand AS seq_CT_strand,
            GA_strand AS seq_GA_strand,
            cigar,
            seq,
            score_before_recal,
            r_strand
          FROM 
            ${DATASET_ID}.${SAMPLE}_recal_sam
          WHERE
            chr = '${CHR}'
          ),
        variants_and_read_id AS (
          SELECT * FROM ${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_chr${CHR}_first
          )
        SELECT 
            snp_id,
            ref, 
            alt,
            pos,
            chr,
            read_start,
            read_end,
            CT_strand,
            GA_strand,
            seq_CT_strand,
            seq_GA_strand,
            r_strand,
            cigar,
            read_id,
            seq,
            score_before_recal 
      FROM sequences
      INNER JOIN variants_and_read_id
      ON read_id_unique = read_id
      "

bq rm -f -t ${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_chr${CHR}_first
