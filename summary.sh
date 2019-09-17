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
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_asm_snp \
    --replace=true \
    "
    SELECT
        chr,
        snp_id,
        dmr_inf,
        dmr_sup,
        ref_reads AS nb_ref_reads,
        alt_reads AS nb_alt_reads,
        effect AS effect_size,
        wilcoxon_pvalue,
        nb_cpg,
        nb_sig_cpg,
        pos_sig_cpg,
        neg_sig_cpg,
        nb_consecutive_asm
    FROM ${DATASET_ID}.${SAMPLE}_dmr_pvalue
    WHERE 
        wilcoxon_pvalue < 0.05 
        AND ABS(effect) > ${DMR_EFFECT}
        AND (pos_sig_cpg >= ${CPG_PER_DMR}
                OR neg_sig_cpg >= ${CPG_PER_DMR})
        AND nb_consecutive_asm > 0
    "

