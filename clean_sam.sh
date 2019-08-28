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
            genome_strand,
            cigar,
            IF(REGEXP_CONTAINS(read_strand, 'CT'), 'READ_1','READ_2') AS insert_read,
            seq
        FROM
            ${DATASET_ID}.${SAMPLE}_recal_sam_raw
        "

