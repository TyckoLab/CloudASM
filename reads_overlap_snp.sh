#!/bin/bash

# Launch the query and save into a temporary table
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_chr${CHR}_${INF}_${SUP} \
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
            score_before_recal,
            r_strand
          FROM 
            ${DATASET_ID}.${SAMPLE}_recal_sam
          WHERE
            chr = '${CHR}'
            AND read_start >= ${INF}
            AND read_end <= ${SUP}
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
            AND pos >= ${INF}
            AND pos <= ${SUP}
      ),
      unique_read_id AS (
        SELECT 
          read_id AS read_id_unique,
          min(read_start) AS seq_start,
          max(read_end) AS seq_end,
          ANY_VALUE(sequences_chr) AS sequences_chr,
          ANY_VALUE(seq_CT_strand) AS seq_CT_strand,
          ANY_VALUE(seq_GA_strand) AS seq_GA_strand
        FROM sequences
        GROUP BY read_id
      ),
      -- Find variants on DNA fragments where R1 or R2 has an overlap
      variants_and_read_id AS (
        SELECT snp_id, pos, ref, alt, read_id_unique, CT_strand, GA_strand 
        FROM variants
        INNER JOIN unique_read_id
        ON 
            seq_start <= pos
            AND seq_end >= pos
            AND (seq_CT_strand = CT_strand OR seq_GA_strand = GA_strand)
      ),
      variants_and_all_reads AS (
      SELECT * FROM variants_and_read_id
      INNER JOIN sequences 
      ON read_id_unique = read_id
      )
      -- Table of all variants x reads where the variant is in the read or its paired read.
      SELECT
        snp_id,
        ref, 
        alt,
        pos,
        sequences_chr AS chr,
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
      FROM variants_and_all_reads
      "
