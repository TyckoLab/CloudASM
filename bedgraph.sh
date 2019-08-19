#!/bin/bash

# Re-order the OB file and transform methylation status in 1 or 0
# Merge the OB and OT
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_both_context \
    --replace=true \
    "WITH OB AS (
            SELECT 
                read_id, 
                chr, 
                pos-1 as pos, 
                IF(meth_call = 'Z',1,0) as meth, 
                1 as cov,
                'OB' as strand
            FROM 
                ${DATASET_ID}.${SAMPLE}_CpGOB
            ),
          OT AS (
            SELECT 
                read_id,
                chr, 
                pos, 
                IF(meth_call = 'Z',1,0) as meth, 
                1 as cov,
                'OT' as strand
            FROM 
                ${DATASET_ID}.${SAMPLE}_CpGOT
            )
          SELECT * FROM OB
          UNION ALL 
          SELECT * FROM OT"

#Sum methylation and coverage per pair of (chr, pos)
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_both_context_bed\
    --replace=true \
    " SELECT 
        chr, 
        pos as pos_start, 
        pos as pos_end, 
        SUM(meth)/SUM(cov) as meth_perc, 
        SUM(cov) as cov
    FROM 
        ${DATASET_ID}.${SAMPLE}_both_context
    GROUP BY 
        chr, pos"

# Methylation percentage bedgraph
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_methyperc_bedgraph \
    --replace=true \
    " SELECT 
        chr, 
        pos_start, 
        pos_end, 
        meth_perc
    FROM 
        ${DATASET_ID}.${SAMPLE}_both_context_bed
    WHERE 
        cov >= 10"

# CpG coverage bedgraph
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_CpGcov_bedgraph \
    --replace=true \
    " SELECT 
        chr, 
        pos_start, 
        pos_end, 
        cov
    FROM 
        ${DATASET_ID}.${SAMPLE}_both_context_bed
    WHERE 
        cov >= 10"

################### Export bedgraph files to bucket
bq extract \
    --field_delimiter "\t" \
    --print_header=false \
    ${DATASET_ID}.${SAMPLE}_methyperc_bedgraph \
    gs://$OUTPUT_B/$SAMPLE/bedgraph/${SAMPLE}_methyperc.bedgraph

bq extract \
    --field_delimiter "\t" \
    --print_header=false \
    ${DATASET_ID}.${SAMPLE}_CpGcov_bedgraph \
    gs://$OUTPUT_B/$SAMPLE/bedgraph/${SAMPLE}_CpGcov.bedgraph


########################## Delete most BQ files

# Delete intermediary sets
#bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_CpGOB
#bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_CpGOT
#bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_tmp

# Delete bedgraph files from BQ
bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_methyperc_bedgraph
bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_CpGcov_bedgraph
