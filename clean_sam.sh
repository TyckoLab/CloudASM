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
            IF(REGEXP_CONTAINS(genome_strand, 'CT'), TRUE, FALSE) AS CT_strand,
            IF(REGEXP_CONTAINS(genome_strand, 'GA'), TRUE, FALSE) AS GA_strand,
            cigar,
            seq,
            score_before_recal
        FROM
            ${DATASET_ID}.${SAMPLE}_recal_sam_raw
        "



## Note the read_strand in the SAM indicates if it is read #1 (CT) or read #2 (GA),
## not to be confused with the genome_strand, which we use in the genotyping.
