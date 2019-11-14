#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_context_filtered \
    --replace=true \
    "
    WITH 
        CONTEXT AS (
            SELECT 
                chr,
                pos,
                meth,
                cov,
                read_id
            FROM ${DATASET_ID}.${SAMPLE}_context
        ),
        SNPS_FOR_CPG AS (
            SELECT 
                chr_snp, 
                pos_snp
            FROM ${DATASET_ID}.${SAMPLE}_snps_for_cpg
        ),
        C_FROM_CONTEXT AS (
            SELECT DISTINCT
                chr, 
                pos AS pos_c 
            FROM CONTEXT
        ),
        REMOVE_C AS (
            SELECT chr, pos_c
            FROM C_FROM_CONTEXT
            EXCEPT DISTINCT 
            SELECT * FROM SNPS_FOR_CPG
        ),
        CHANGE_TO_G AS (
            SELECT 
                chr AS chr_clean, 
                pos_c + 1 AS pos_g 
            FROM REMOVE_C
        ),
        REMOVE_C_AND_G AS (
        SELECT * FROM CHANGE_TO_G
        EXCEPT DISTINCT
        SELECT * FROM SNPS_FOR_CPG
        )
        SELECT chr, pos, meth, cov, read_id 
        FROM CONTEXT INNER JOIN REMOVE_C_AND_G
        ON chr = chr_clean AND pos = pos_g - 1
        "

# Create a long table of all CpG x SNP x read_id combinations
# 90% of CpGs are dropped because their respective reads could not be linked to a REF or ALT of any SNP
# This table will be used in the DMR calculation.
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_cpg_read_genotype \
    --replace=true \
    "
    WITH 
        CONTEXT AS (
            SELECT * 
            FROM ${DATASET_ID}.${SAMPLE}_context_filtered
        ),
        GENOTYPE AS (
            SELECT 
                snp_id,
                pos AS snp_pos,
                read_id AS geno_read_id,
                allele
            FROM ${DATASET_ID}.${SAMPLE}_vcf_reads_genotype
        ),
        -- create a table with each combination of snp, CpG, read_id, REF, ALT
        COMBINED AS (
            SELECT * FROM CONTEXT
            INNER JOIN GENOTYPE 
            ON read_id = geno_read_id
        )
        -- we remove the extra columns and keep distinct rows
        SELECT
            chr, 
            pos,
            meth,
            cov,
            snp_id,
            allele,
            read_id
        FROM COMBINED
        "

# This file will be used to compute single CPG ASM
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_cpg_genotype \
    --replace=true \
        "
        WITH 
        ----------------------------------------------------------
        -- Find CpGs that are at least covered 5x on both REF and ALT for a SNP
        ----------------------------------------------------------
        CLEAN AS (
            SELECT * FROM ${DATASET_ID}.${SAMPLE}_cpg_read_genotype
        ),
        REF_ALLELES AS (
            SELECT 
                chr,
                pos,
                snp_id,
                SUM(meth) AS ref_meth,
                SUM(cov) AS ref_cov
            FROM CLEAN
            WHERE allele = 'REF'
            GROUP BY chr, pos, snp_id
        ),
        ALT_ALLELES AS (
            SELECT 
                chr AS chr_alt,
                pos AS pos_alt,
                snp_id AS snp_id_alt,
                SUM(meth) AS alt_meth,
                SUM(cov) AS alt_cov
            FROM CLEAN
            WHERE allele = 'ALT'
            GROUP BY chr_alt, pos_alt, snp_id_alt
        ),
        REF_AND_ALT AS (
            SELECT * FROM REF_ALLELES
            INNER JOIN ALT_ALLELES
            ON chr = chr_alt AND 
                pos = pos_alt AND
                snp_id = snp_id_alt
        )
            SELECT 
                chr,
                pos,
                snp_id,
                ref_cov,
                ref_meth,
                alt_cov,
                alt_meth
            FROM REF_AND_ALT
            -- we require that each allele is covered 5x
            WHERE ref_cov >=${CPG_COV} AND alt_cov >= ${CPG_COV}
        "

# Save file in bucket to use it in a Python script to compute a Fisher exact test
bq extract \
    --field_delimiter "," \
    --print_header=true \
    ${DATASET_ID}.${SAMPLE}_cpg_genotype \
    gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_cpg_genotype.csv