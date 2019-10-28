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
          r_strand,
          cigar,
          read_id,
          seq,
          score_before_recal
        FROM VARIANTS_SEQUENCES
      ), 
      -- Now we need to recover the reads that do not overlap the variant but sequence the same DNA fragment.
      -- First identify DNA segments where both reads overlap the variant position
      READ_ID_WITH_SNP_ON_BOTH_STRANDS AS (
        SELECT COUNT(r_strand), read_id AS read_id_both_strands
        FROM VARIANTS_SEQUENCES
        GROUP BY read_id   
        HAVING COUNT(r_strand) > 1
      ),
      -- Then identify (read_id, strand) that do not overlap the variant position but are on 
      -- a DNA segment overlapping the variant
      READ_STRANDS_WITH_NO_SNP_BUT_OVERLAPPING_DNA_SEGMENT_WITH_SNP AS (
        SELECT 
          snp_id, 
          FALSE AS snp_in_read,
          chr,
          pos,
          ref,
          alt,
          CT_strand,
          GA_strand,
          read_id AS read_id_to_find, 
          IF(REGEXP_CONTAINS(r_strand, 'R1'), 'R2', 'R1') AS r_strand_new 
        FROM VARIANTS_SEQUENCES
        FULL JOIN READ_ID_WITH_SNP_ON_BOTH_STRANDS
        ON read_id = read_id_both_strands
        WHERE read_id_both_strands IS NULL
      ),
      SEQUENCES_WITHOUT_SNP_BUT_OVERLAPPING_DNA_SEGMENT_WITH_SNP AS (
        SELECT * FROM READ_STRANDS_WITH_NO_SNP_BUT_OVERLAPPING_DNA_SEGMENT_WITH_SNP
        INNER JOIN sequences
        ON read_id_to_find = read_id AND r_strand_new = r_strand
      )
      SELECT 
        snp_id,
        snp_in_read,
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
        r_strand,
        cigar,
        read_id,
        seq,
        score_before_recal
       FROM VARIANTS_AND_OVERLAPPING_READS
      UNION ALL
      SELECT 
          snp_id,
          snp_in_read,
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
          r_strand,
          cigar,
          read_id,
          seq,
          score_before_recal
        FROM SEQUENCES_WITHOUT_SNP_BUT_OVERLAPPING_DNA_SEGMENT_WITH_SNP
      "



