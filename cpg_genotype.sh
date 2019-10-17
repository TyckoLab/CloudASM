#!/bin/bash

# We remove all CpGs where the C or G overlap with a SNP found in the "raw" Bis-SNP database
# 60 seconds for 2 chromosomes
# That excludes about 5% of all CGs.

bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_context_filtered \
    --replace=true \
    "
    WITH 
        CONTEXT AS (
            SELECT * 
            FROM ${DATASET_ID}.${SAMPLE}_context
        ),
        VCF_RAW AS (
            SELECT 
                chr_snp, 
                pos_snp AS pos_snp_c, 
                pos_snp + 1 AS pos_snp_g
            FROM ${DATASET_ID}.${SAMPLE}_snps_for_cpg
        ),
        CONTEXT_FILTERED_C AS (
            SELECT distinct chr, pos FROM CONTEXT
            EXCEPT DISTINCT
            SELECT distinct chr_snp, pos_snp_c FROM VCF_RAW
        ),
        CONTEXT_FILTERED_CG AS (
            SELECT chr, pos FROM CONTEXT_FILTERED_C
            EXCEPT DISTINCT
            SELECT chr_snp, pos_snp_g FROM VCF_RAW
        ),
        FILTERED AS (
            SELECT chr AS chr_filt, pos AS pos_filt
            FROM CONTEXT_FILTERED_CG
        ),
        CONTEXT_FILTERED_READ AS (
        SELECT * FROM FILTERED
        INNER JOIN CONTEXT
        ON chr_filt = chr AND pos = pos_filt
        )
        SELECT chr, pos, meth, cov, read_id FROM CONTEXT_FILTERED_READ
        --SELECT * FROM VCF_RAW
        "

# Create a table with as many rows as possible for CpG x SNP x read_id
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
        SELECT DISTINCT
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