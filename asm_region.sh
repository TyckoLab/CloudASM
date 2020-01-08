#!/bin/bash

# Import the file where the p-value of each cpg
# is calculated.
bq --location=US load \
               --replace=true \
               --source_format=CSV \
               --field_delimiter "," \
                ${DATASET_ID}.${SAMPLE}_cpg_asm \
               gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_cpg_asm.csv \
               chr:STRING,pos:INTEGER,snp_id:STRING,snp_pos:INTEGER,ref_cov:INTEGER,ref_meth:INTEGER,alt_cov:INTEGER,alt_meth:INTEGER,fisher_pvalue:FLOAT


# Delete SAMPLE_cpg_genotype (the new file sample_cpg_asm has the same information and the ASM p-value)
bq rm -f -t ${DATASET_ID}.${SAMPLE}_cpg_genotype

# Create a table with one row per SNP with at least CPG_PER_ASM_REGION CpGs nearby for which a fisher p-value was calculated
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_snp_cpg_array \
    --replace=true \
    "
    WITH 
    -- SNPs with their respective arrays of CpGs (each coming with their fisher p-value)
        SNP_CPG_ARRAY AS (
            SELECT
                snp_id,
                ANY_VALUE(snp_pos) AS snp_pos,
                ANY_VALUE(chr) AS chr,
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
                snp_id AS snp_id_dmr, 
                snp_pos,
                chr AS chr_dmr, 
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
        )
        -- Keep DMR with at least CPG_PER_ASM_REGION CpGs. For 3 CpG per DMR, half of the SNPs are dropped
            SELECT * FROM HET_SNP WHERE nb_cpg >= ${CPG_PER_ASM_REGION}
        "

# Create a table of all CpGs to be used in DMR effect
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_cpg_for_dmr_effect \
    --replace=true \
    "
    WITH 
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
        )
        -- Keep the combination of CpG, read_id, allele for which CpGs are at least 5x covered on each allele
            SELECT DISTINCT
                chr_cpg,
                pos_cpg,
                meth,
                cov,
                snp_id_cpg,
                allele,
                read_id  
            FROM ALL_CPG
            INNER JOIN WELL_COVERED_CPG
            ON chr_cpg = chr_well_cov 
               AND pos_cpg = pos_well_cov 
        "

# Create a table of SNPs and an array of reads and their fractional methylation 
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_cpg_for_wilcox \
    --replace=true \
    "
    WITH        
        HET_SNP AS (
            SELECT * 
            FROM ${DATASET_ID}.${SAMPLE}_snp_cpg_array
        ),
        CPG_FOR_DMR AS (
            SELECT *
            FROM ${DATASET_ID}.${SAMPLE}_cpg_for_dmr_effect
        )
        -- we take all CpGs across the DMR that are found in between 2 significant ASM CpGs if they exist, 
        -- anywhere in the DMR otherwise (which will not be a real DMR because we demand significant CpGs in a DMR).
        -- We do not discard them to compute accurately the corrected Wilcoxon p-value.
        SELECT
            chr_cpg, 
            pos_cpg, 
            meth, 
            cov, 
            snp_id_cpg AS snp_id, 
            snp_pos,
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
        FROM HET_SNP
        INNER JOIN CPG_FOR_DMR
        ON 
            snp_id_cpg = snp_id_dmr 
            AND chr_cpg = chr_dmr
        WHERE 
            pos_cpg >= IF(dmr_inf IS NULL, min_cpg, dmr_inf) 
            AND pos_cpg <= IF(dmr_sup IS NULL, max_cpg, dmr_sup)
        "
        
bq query \
    --use_legacy_sql=false \
    --destination_table ${DATASET_ID}.${SAMPLE}_snp_for_dmr \
    --replace=true \
    "
    WITH 
        QUALIFYING_CPG_WILCOX AS (
            SELECT * FROM ${DATASET_ID}.${SAMPLE}_cpg_for_wilcox
        ),
        METHYL_PER_READ_WILCOX AS (
            SELECT 
                snp_id,
                snp_pos,
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
            GROUP BY snp_id, snp_pos, read_id, allele, chr_cpg
        ),
        SNP_METHYL_ARRAY_REF_WILCOX AS (
            SELECT
                snp_id,
                ANY_VALUE(snp_pos) AS snp_pos,
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
        SELECT 
            snp_id,
            snp_pos,
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


# Delete intermediary files
bq rm -f -t ${DATASET_ID}.${SAMPLE}_snp_cpg_array
bq rm -f -t ${DATASET_ID}.${SAMPLE}_cpg_for_dmr_effect
bq rm -f -t ${DATASET_ID}.${SAMPLE}_cpg_for_wilcox

# Export file to JSON format in the bucket
# (nested arrays are not supported in)
bq extract \
    --destination_format NEWLINE_DELIMITED_JSON \
    ${DATASET_ID}.${SAMPLE}_snp_for_dmr \
    gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_snp_for_dmr.json