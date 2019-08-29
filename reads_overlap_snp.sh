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
            insert_length,
            seq_length,
            cigar,
            insert_read,
            seq
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
          cov,
          CT_strand,
          GA_strand
        FROM 
          ${DATASET_ID}.${SAMPLE}_vcf
        WHERE 
            chr = '${CHR}'
      )
      SELECT *
      FROM 
        variants
      INNER JOIN
        sequences ON
          sequences.read_start <= pos
          AND sequences.read_end >= pos
          AND sequences.sequences_chr = chr
          AND (sequences.seq_CT_strand = CT_strand OR sequences.seq_GA_strand = GA_strand)
    "

# # Save the table in Google Cloud Storage

# bq extract \
#     --field_delimiter "," \
#     --print_header=true \
#     ${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_${CHR} \
#     gs://$OUTPUT_B/$SAMPLE/tmp/${SAMPLE}_vcf_reads_tmp_${CHR}.vcf
