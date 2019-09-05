#!/bin/bash

# Clean the SAM
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_recal_sam \
    --replace=true \
        "SELECT
            read_id,
            chr,
            read_start,
            read_start + BYTE_LENGTH(seq) -1 AS read_end,
            length AS insert_length,
            BYTE_LENGTH(seq) AS seq_length,
            -- CT_strand is read 1 and GA_strand is read 2
            IF(REGEXP_CONTAINS(genome_strand, 'CT'), TRUE, FALSE) AS CT_strand,
            IF(REGEXP_CONTAINS(genome_strand, 'GA'), TRUE, FALSE) AS GA_strand,
            cigar,
            seq,
            score_before_recal
        FROM
            ${DATASET_ID}.${SAMPLE}_recal_sam_raw
        "

