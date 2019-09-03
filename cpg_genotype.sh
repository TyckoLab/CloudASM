#!/bin/bash

bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_cpg_genotype \
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
                read_id AS geno_read_id,
                allele 
            FROM ${DATASET_ID}.${SAMPLE}_vcf_reads_genotype
        ),
        -- create a table with each combination of snp, CpG, REF, ALT
        COMBINED AS (
            SELECT * FROM CONTEXT
            INNER JOIN GENOTYPE 
            ON read_id = geno_read_id
        ),
        -- we remove the extra columns in CLEAN
        CLEAN AS (
            SELECT 
                chr, 
                pos,
                meth,
                cov,
                snp_id,
                allele,
                read_id
            FROM COMBINED
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
       -- ref_meth - alt_meth AS m_ref_minus_alt,
       -- ROUND(SAFE_DIVIDE(ref_meth - alt_meth,ref_meth+alt_meth),3) AS meth_perc
    FROM REF_AND_ALT
    -- we require that each allele is covered 5x
    WHERE ref_cov >=5 AND alt_cov >= 5
    "

# Save file in bucket to use it in a Python script to compute a Fisher exact test

bq extract \
    --field_delimiter "," \
    --print_header=true \
    ${DATASET_ID}.${SAMPLE}_cpg_genotype \
    gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_cpg_genotype.csv