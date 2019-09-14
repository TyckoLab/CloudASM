#!/bin/bash

# This file will be used to compute the DMRs
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_cpg_read_genotype \
    --replace=true \
    "
    WITH 
        CONTEXT AS (
            SELECT * 
            FROM ${DATASET_ID}.${SAMPLE}_context
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
        -- we remove the extra columns
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
        ),
        REF_AND_ALT_CLEAN AS (
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
        ),
        ----------------------------------------------------------
        -- keep CpG that where the C does not overlap with any SNP
        ----------------------------------------------------------
        CPG_NO_SNP_ON_C AS (
            SELECT chr,pos 
                FROM REF_AND_ALT_CLEAN
            EXCEPT DISTINCT
            SELECT chr,pos 
                FROM ${DATASET_ID}.${SAMPLE}_vcf_reads_genotype
        ),
        CPG_NO_SNP_ON_C_AND_G AS (
            SELECT chr,pos 
                FROM CPG_NO_SNP_ON_C
            EXCEPT DISTINCT
            SELECT chr, pos - 1 AS pos 
                FROM ${DATASET_ID}.${SAMPLE}_vcf_reads_genotype
        ),
        REF_AND_ALT_CLEAN_FORMATED AS (
            SELECT
                chr AS chr_tmp,
                pos AS pos_tmp,
                snp_id,
                ref_cov,
                ref_meth,
                alt_cov,
                alt_meth
            FROM REF_AND_ALT_CLEAN 
        ),
        CPG_GENOTYPE AS (
            SELECT * 
                FROM CPG_NO_SNP_ON_C_AND_G 
            INNER JOIN 
            REF_AND_ALT_CLEAN_FORMATED
            ON chr = chr_tmp AND pos = pos_tmp
        )
        SELECT
            chr,
            pos,
            snp_id,
            ref_cov,
            ref_meth,
            alt_cov,
            alt_meth
            FROM CPG_GENOTYPE
    "

# Save file in bucket to use it in a Python script to compute a Fisher exact test
bq extract \
    --field_delimiter "," \
    --print_header=true \
    ${DATASET_ID}.${SAMPLE}_cpg_genotype \
    gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_cpg_genotype.csv