#!/bin/bash

# This process will filter out about 30% of all CpG 
#because they do not have at least 20% methylation or 10x coverage

# Re-order the OB file and transform methylation status in 1 or 0
echo "Merge the OB and OT"
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_both_context_tmp \
    --replace=true \
    "WITH 
        OB AS (
            SELECT 
                read_id, 
                chr, 
                -- Bismark gives the position of the C in the negative strand. need to 
                -- subtract 1 to have the position of the CG in positive strand.
                pos-1 as pos, 
                IF(meth_call = 'Z',1,0) as meth, 
                1 as cov
            FROM 
                ${DATASET_ID}.${SAMPLE}_CpGOB
            ),
        OT AS (
            SELECT 
                read_id,
                chr, 
                pos, 
                IF(meth_call = 'Z',1,0) as meth, 
                1 as cov
            FROM 
                ${DATASET_ID}.${SAMPLE}_CpGOT
            )
          SELECT * FROM OB
          UNION ALL 
          SELECT * FROM OT"

echo "Sum methylation and coverage per CpG. Keep at least 10x cov."
# We impose a coverage of at least 10x per CpG.
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_merged_context_bed \
    --replace=true \
    "  
    WITH MERGED_STRANDS AS (
        SELECT 
            chr, 
            pos AS pos_start, 
            pos AS pos_end, 
            SUM(cov) AS cov,
            SUM(meth)/SUM(cov) AS meth_perc
        FROM 
            ${DATASET_ID}.${SAMPLE}_both_context_tmp
        GROUP BY 
            chr, pos
        )
        SELECT 
            chr,
            pos_start,
            pos_end,
            cov,
            meth_perc
        FROM 
            MERGED_STRANDS
        WHERE
            -- 2 x the minimum coverage required by allele (no allele distinction at this stage)
            cov >= 2*${CPG_COV}
    "

echo "Filter out from both_context_tmp the CpG sites that do not have at least 10x cov"
# This removes about 20% of all CpG flagged by Bismark.
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_context \
    --replace=true \
    " 
    WITH 
        MERGED_CONTEXT_FILT AS (
            SELECT
                chr AS chr_cpg,
                pos_start AS pos_cpg
            FROM
                ${DATASET_ID}.${SAMPLE}_merged_context_bed
        )
        SELECT 
            read_id,
            chr,
            pos,
            meth,
            cov
        FROM
            ${DATASET_ID}.${SAMPLE}_both_context_tmp
        INNER JOIN
            MERGED_CONTEXT_FILT ON
            chr_cpg = chr
            AND pos_cpg = pos
    "

echo "Methylation percentage bedgraph"
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
        ${DATASET_ID}.${SAMPLE}_merged_context_bed
    "

echo "CpG coverage bedgraph"
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
        ${DATASET_ID}.${SAMPLE}_merged_context_bed
    "

################### Export bedgraph files to bucket

echo "Export bedgraph files of net methylation and coverage"

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

echo "Delete bedgraph files from BQ"
bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_both_context_tmp
bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_methyperc_bedgraph
bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_CpGcov_bedgraph
bq rm -f -t ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_merged_context_bed

echo "Delete raw files"
bq rm -f -t ${DATASET_ID}.${SAMPLE}_CpGOB
bq rm -f -t ${DATASET_ID}.${SAMPLE}_CpGOT