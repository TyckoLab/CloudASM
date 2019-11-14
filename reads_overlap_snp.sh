#!/bin/bash

# Build a table where each combination of read and variant where the variant is overlapping the read
# This will discard all paired reads not overlapping the variant but for which we know the genotyping
# From the first read. The inner join will be performed in the next script
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_chr${CHR} \
    --replace=true \
    " WITH
        sequences AS (
          SELECT 
            read_id,
            read_start,
            read_end,
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
      )
      -- Find variants on DNA fragments where R1 or R2 has an overlap
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
          r_strand,
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
