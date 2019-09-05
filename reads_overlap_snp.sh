#!/bin/bash

# Launch the query and save into a temporary table
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_${CHR} \
    --replace=true \
    " WITH
        sequences AS (
          SELECT 
            read_id,
            read_start,
            read_end,
            chr AS sequences_chr,
            CT_strand AS seq_CT_strand,
            GA_strand AS seq_GA_strand,
            cigar,
            seq,
            score_before_recal
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
      VARIANTS_SEQUENCES AS (
        SELECT * FROM variants
        INNER JOIN
          sequences ON
            sequences.read_start <= pos
            AND sequences.read_end >= pos
            AND sequences.sequences_chr = chr
            AND (sequences.seq_CT_strand = CT_strand OR sequences.seq_GA_strand = GA_strand)
      )
      SELECT 
        snp_id,
        chr,
        pos,
        ref,
        alt,
        read_start,
        read_end,
        CT_strand,
        GA_strand,
        seq_CT_strand,
        seq_GA_strand,
        cigar,
        read_id,
        seq,
        score_before_recal
      FROM VARIANTS_SEQUENCES
    "
