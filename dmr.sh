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


# Delete SAMPLE_cpg_genotype (the new file sample_cpg_asm has the same information and the ASM p-value)
bq rm -f -t ${DATASET_ID}.${SAMPLE}_cpg_genotype

# Query to select the SNPs with at least 3 significant CpGs in the same direction
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_snp_for_dmr\
    --replace=true \
    "
    WITH 
        DMR_BOUNDS AS (
            SELECT
                snp_id,
                ANY_VALUE(chr) AS chr,
                --ANY_VALUE(nb_cpg) AS nb_cpg,
                ARRAY_AGG(STRUCT(pos,
                                ROUND(alt_meth/alt_cov-ref_meth/ref_cov,3) AS alt_minus_ref,
                                fisher_pvalue)) 
                        AS cpg
            FROM ${DATASET_ID}.${SAMPLE}_cpg_asm
            GROUP BY snp_id
        ),
        DMR_LIST AS (
            SELECT snp_id, 
                chr, 
                ARRAY_LENGTH(cpg) AS nb_cpg, 
                (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg) WHERE fisher_pvalue < 0.05) AS nb_sig_cpg, 
                (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg) WHERE fisher_pvalue < 0.05 AND SIGN(alt_minus_ref) = 1) AS pos_sig_cpg,
                (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg) WHERE fisher_pvalue < 0.05 AND SIGN(alt_minus_ref) = -1) AS neg_sig_cpg, 
                (SELECT MIN(pos) FROM UNNEST(cpg) WHERE fisher_pvalue < 0.05) AS min_cpg,
                (SELECT MAX(pos) FROM UNNEST(cpg) WHERE fisher_pvalue < 0.05) AS max_cpg,
                cpg
            FROM DMR_BOUNDS
        ),
        -- INTERSECT WITH THE LIST OF CPGS WHERE WE HAVE THE SNP AND READ_ID, AND ALLELE
       -- SELECT snp_id AS snp_id_for_dmr, chr, nb_cpg, min_cpg, max_cpg FROM DMR_LIST WHERE pos_sig_cpg >=3 OR neg_sig_cpg >=3
        SELECT * FROM DMR_LIST
    "

# Intersect this list of SNPs with the the CpG files (sample_cpg_read_genotype)
# The goal is to create a table where each SNP comes with an array
