#!/bin/bash

# Import the file where the p-value of each cpg
# is calculated.
bq --location=US load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "," \
                ${DATASET_ID}.${SAMPLE}_cpg_asm \
               gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_cpg_asm.csv \
               chr:STRING,pos:INTEGER,snp_id:STRING,ref_cov:INTEGER,ref_meth:INTEGER,alt_cov:INTEGER,alt_meth:INTEGER,fisher_pvalue:FLOAT


bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_tmp\
    --replace=true \
    "WITH 
        NB_CPG_PER_SNP AS (
            SELECT snp_id, COUNT(*) AS nb_cpg
            FROM ${DATASET_ID}.${SAMPLE}_cpg_asm
            GROUP BY snp_id 
            ),
        -- all SNPS with at least 3 CpGs with ASM. They may not be significant
        SNP_WITH_POTENTIAL_DMR AS (
            SELECT * FROM NB_CPG_PER_SNP
            WHERE nb_cpg >=3
        )
        SELECT * FROM SNP_WITH_POTENTIAL_DMR
    "

# WITH TEST AS (
#   SELECT  *  FROM `hackensack-tyco.wgbs_asm.gm12878_cpg_asm` 
#   )
#  SELECT TO_JSON_STRING(TEST) FROM TEST