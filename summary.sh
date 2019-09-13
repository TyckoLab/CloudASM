#!/bin/bash

# Import the file where the p-value of each cpg
# is calculated.
bq --location=US load \
               --autodetect \
               --replace=true \
               --source_format=NEWLINE_DELIMITED_JSON \
                ${DATASET_ID}.${SAMPLE}_dmr_pvalue \
               gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_dmr_pvalue.json 
               


# Delete the file generated before computing the p-values.
bq rm -f -t ${DATASET_ID}.${SAMPLE}_snp_for_dmr

# Query to select the SNPs with at least 3 significant CpGs in the same direction
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.asm_snp \
    --replace=false \
    "
    SELECT
        '${SAMPLE}' as sample,
        chr,
        snp_id,
        ref_reads AS nb_ref_reads,
        alt_reads AS nb_alt_reads,
        effect AS effect_size,
        wilcoxon_pvalue
    FROM ${DATASET_ID}.${SAMPLE}_dmr_pvalue
    WHERE 
        wilcoxon_pvalue < 0.05
 
    "

