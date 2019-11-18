#!/bin/bash

# Build a table where each combination of read and variant where the variant is overlapping the read
# This will discard all paired reads not overlapping the variant but for which we know the genotyping

# Extract from the SAM file
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_vcf_reads_chr${CHR} \
    --replace=true \
    "
     WITH 
      full_sam AS (
        SELECT 
          read_id, 
          chr, 
          read_start, 
          read_end, 
          CT_strand AS seq_CT_strand, 
          GA_strand AS seq_GA_strand, 
          r_strand, 
          cigar, 
          seq, 
          score_before_recal 
        FROM ${DATASET_ID}.${SAMPLE}_recal_sam
        WHERE chr = '${CHR}'
      ),
      trimmed_sam AS (
        SELECT 
          read_start AS t_start,
          read_end AS t_end,
          read_id AS t_read,
          r_strand AS t_strand
        FROM full_sam
      ),
      vcf AS (
        SELECT 
          snp_id, 
          pos, 
          ref, 
          alt, 
          CT_strand, 
          GA_strand 
        FROM ${DATASET_ID}.${SAMPLE}_vcf
        WHERE chr = '${CHR}'
      ),
      -- Perform a Broadcast join with the minimum amount of columns
      sam_vcf_join AS (
        SELECT * FROM trimmed_sam 
        INNER JOIN vcf 
        ON (pos >= t_start AND pos <= t_end)
      )
      -- Perform a Hash join to get the other required columns.
      SELECT
        chr,
        snp_id,
        pos, ref,
        alt,
        CT_strand,
        GA_strand,
        seq_CT_strand,
        seq_GA_strand,
        read_start,
        read_end,
        cigar,
        read_id,
        seq,
        score_before_recal
      FROM full_sam 
      INNER JOIN sam_vcf_join 
      ON t_read = read_id AND r_strand = t_strand
      "
