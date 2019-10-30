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


echo "Make a table with all the possible DMRs (before computing p-value)"
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_snp_for_dmr \
    --replace=true \
    "
    WITH 
    -- SNPs with their respective arrays of CpGs (each coming with their fisher p-value)
        SNP_CPG_ARRAY AS (
            SELECT
                snp_id,
                ANY_VALUE(chr) AS chr,
                --ANY_VALUE(nb_cpg) AS nb_cpg,
                ARRAY_AGG(
                    STRUCT(
                        pos,
                        ROUND(alt_meth/alt_cov-ref_meth/ref_cov,3) AS effect,
                        fisher_pvalue,
                        ref_cov,
                        alt_cov
                    )
                    ORDER BY pos
                ) AS cpg
            FROM ${DATASET_ID}.${SAMPLE}_cpg_asm
            GROUP BY snp_id
        ),
        -- Extract the number of CpG per SNP 
        HET_SNP AS (
            SELECT 
                snp_id, 
                chr, 
                ARRAY_LENGTH(cpg) AS nb_cpg, 
                (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg) WHERE fisher_pvalue < ${P_VALUE}) AS nb_sig_cpg, 
                (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg) WHERE fisher_pvalue < ${P_VALUE} AND SIGN(effect) = 1) AS pos_sig_cpg,
                (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg) WHERE fisher_pvalue < ${P_VALUE} AND SIGN(effect) = -1) AS neg_sig_cpg, 
                (SELECT MIN(pos) FROM UNNEST(cpg) WHERE fisher_pvalue < ${P_VALUE}) AS dmr_inf,
                (SELECT MAX(pos) FROM UNNEST(cpg) WHERE fisher_pvalue < ${P_VALUE}) AS dmr_sup,
                (SELECT MIN(pos) FROM UNNEST(cpg)) AS min_cpg,
                (SELECT MAX(pos) FROM UNNEST(cpg)) AS max_cpg,
                cpg
            FROM SNP_CPG_ARRAY
        ),
        -- For 3 CpG per DMR, half of the SNPs are dropped
        DMR_WITH_ENOUGH_CPG AS (
            SELECT * FROM HET_SNP WHERE nb_cpg >= ${CPG_PER_DMR}
        ),  
        SNP_FOR_DMR AS (
            SELECT
                snp_id AS snp_id_dmr, 
                chr AS chr_dmr,
                min_cpg,
                max_cpg,
                dmr_inf,
                dmr_sup,
                nb_cpg,
                nb_sig_cpg,
                pos_sig_cpg,
                neg_sig_cpg,
                cpg
            FROM DMR_WITH_ENOUGH_CPG
        ),
        -- Import the list of CpGs with their genotype as a function of read_id and snp_id
        ALL_CPG AS (
            SELECT 
                chr AS chr_cpg,
                pos AS pos_cpg,
                meth,
                cov,
                snp_id AS snp_id_cpg,
                allele,
                read_id 
            FROM ${DATASET_ID}.${SAMPLE}_cpg_read_genotype
        ),
        -- Import all CpGs that are at least 5x covered on each allele and for which we have a fisher p-value
        WELL_COVERED_CPG AS (
            SELECT 
                chr AS chr_well_cov, 
                pos AS pos_well_cov
            FROM ${DATASET_ID}.${SAMPLE}_cpg_asm
        ),
        -- Keep the combination of CpG, read_id, allele for which CpGs are at least 5x covered on each allele
        FILTERED_CPG AS (
            SELECT * FROM ALL_CPG
            INNER JOIN WELL_COVERED_CPG
            ON chr_cpg = chr_well_cov 
               AND pos_cpg = pos_well_cov 
        ),
        -- same, but keep distinct rows
        FILTERED_CPG_CLEAN AS (
            SELECT DISTINCT
                chr_cpg,
                pos_cpg,
                meth,
                cov,
                snp_id_cpg,
                allele,
                read_id 
            FROM FILTERED_CPG
        ),
        -- Creates a long table of snps (with their arrays of CpGs) and CpG (with their read_id/genotype) 
        CPG_DMR AS (
            SELECT * FROM SNP_FOR_DMR
            INNER JOIN FILTERED_CPG_CLEAN
            ON snp_id_cpg = snp_id_dmr AND chr_cpg = chr_dmr
        ),
        -- we take all CpGs across the DMR that are found in between 2 significant ASM CpGs if they exist, anywhere in the DMR otherwise (which will not be a real DMR because we demand significant CpGs in a DMR).
        -- We do not discard them to compute accurately the corrected Wilcoxon p-value.
        QUALIFYING_CPG_WILCOX AS (
        -- All CpGs that are located within the boundaries of the potential DMR
        -- This removes 1/3 of all CpGs.
            SELECT 
                chr_cpg, 
                pos_cpg, 
                meth, 
                cov, 
                snp_id_cpg AS snp_id, 
                min_cpg,
                max_cpg,
                dmr_inf,
                dmr_sup,
                allele,
                read_id,
                nb_cpg,
                nb_sig_cpg,
                pos_sig_cpg,
                neg_sig_cpg,
                cpg
            FROM CPG_DMR
            -- Many snps from HET SNPs do not have lower or upper bound.
            WHERE 
                pos_cpg >= IF(dmr_inf IS NULL, min_cpg, dmr_inf) 
                AND pos_cpg <= IF(dmr_sup IS NULL, max_cpg, dmr_sup)
        ),
        METHYL_PER_READ_WILCOX AS (
            SELECT 
                snp_id,
                chr_cpg,
                ANY_VALUE(dmr_inf) AS dmr_inf,
                ANY_VALUE(dmr_sup) AS dmr_sup,
                read_id,
                allele,
                ROUND(SAFE_DIVIDE(SUM(meth),SUM(cov)),5) AS methyl_perc,
                ANY_VALUE(nb_cpg) AS nb_cpg,
                ANY_VALUE(nb_sig_cpg) AS nb_sig_cpg,
                ANY_VALUE(pos_sig_cpg) AS pos_sig_cpg,
                ANY_VALUE(neg_sig_cpg) AS neg_sig_cpg,
                ANY_VALUE(cpg) AS cpg
            FROM QUALIFYING_CPG_WILCOX
            GROUP BY snp_id, read_id, allele, chr_cpg
        ),
        SNP_METHYL_ARRAY_REF_WILCOX AS (
            SELECT
                snp_id,
                ANY_VALUE(chr_cpg) AS chr,
                ANY_VALUE(dmr_inf) AS dmr_inf,
                ANY_VALUE(dmr_sup) AS dmr_sup,
                ARRAY_AGG(STRUCT(methyl_perc)) AS ref,
                ANY_VALUE(nb_cpg) AS nb_cpg,
                ANY_VALUE(nb_sig_cpg) AS nb_sig_cpg,
                ANY_VALUE(pos_sig_cpg) AS pos_sig_cpg,
                ANY_VALUE(neg_sig_cpg) AS neg_sig_cpg,
                ANY_VALUE(cpg) AS cpg
            FROM METHYL_PER_READ_WILCOX
            WHERE allele = 'REF'
            GROUP BY snp_id
        ),
        SNP_METHYL_ARRAY_ALT_WILCOX AS (
            SELECT
                snp_id AS snp_id_alt,
                ARRAY_AGG(STRUCT(methyl_perc)) AS alt
            FROM METHYL_PER_READ_WILCOX
            WHERE allele = 'ALT'
            GROUP BY snp_id
        ),
        SNP_METHYL_JOIN_WILCOX AS (
            SELECT * FROM SNP_METHYL_ARRAY_REF_WILCOX
            INNER JOIN SNP_METHYL_ARRAY_ALT_WILCOX
            ON snp_id = snp_id_alt
        ),
        SNP_METHYL_EFFECT AS (
            SELECT
                snp_id AS snp_id_effect,
                ROUND(((SELECT AVG(methyl_perc) FROM UNNEST(alt)) - (SELECT AVG(methyl_perc) FROM UNNEST(ref))),3) AS effect
            FROM SNP_METHYL_JOIN_WILCOX
        )
        SELECT snp_id,   
               chr,
               dmr_inf,
               dmr_sup,
               effect AS dmr_effect,
               ARRAY_LENGTH(ref) AS ref_reads, 
               ARRAY_LENGTH(alt) AS alt_reads,
               ref, 
               alt,
               nb_cpg,
               nb_sig_cpg,
               pos_sig_cpg,
               neg_sig_cpg,
               cpg 
        FROM SNP_METHYL_EFFECT
        INNER JOIN SNP_METHYL_JOIN_WILCOX
        ON snp_id = snp_id_effect
        "


# Export file to JSON format in the bucket
# (nested arrays are not supported in)
bq extract \
    --destination_format NEWLINE_DELIMITED_JSON \
    ${DATASET_ID}.${SAMPLE}_snp_for_dmr \
    gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_snp_for_dmr.json