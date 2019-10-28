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
            score_before_recal
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
      VARIANTS_SEQUENCES AS (
        SELECT * FROM variants INNER JOIN sequences 
        ON
            sequences.read_start <= pos
            AND sequences.read_end >= pos
            AND sequences.sequences_chr = chr
            AND (sequences.seq_CT_strand = CT_strand OR sequences.seq_GA_strand = GA_strand)
      ),
      VARIANTS_AND_OVERLAPPING_READS AS (
        SELECT 
          snp_id,
          TRUE AS snp_in_read,
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
      ), 
      -- Now we need to recover the reads that do not overlap the variant but sequence the same DNA fragment.
      READ_ID_OVERLAPPING_VARIANTS AS (
        SELECT 
          snp_id,
          read_id AS read_id_to_keep,
          FALSE AS snp_id_read
        FROM VARIANTS_AND_OVERLAPPING_READS
      ),
      ALL_READS AS (
        SELECT * FROM READ_ID_OVERLAPPING_VARIANTS INNER JOIN sequences
        ON read_id_to_keep = read_id AND 
      )
      SELECT * FROM READ_ID_OVERLAPPING_VARIANTS
    "
