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

echo "Query to select the SNPs with at least CPG_PER_DMR nearby."
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_het_snp \
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
                (SELECT MIN(pos) FROM UNNEST(cpg) WHERE fisher_pvalue < ${P_VALUE}) AS min_cpg,
                (SELECT MAX(pos) FROM UNNEST(cpg) WHERE fisher_pvalue < ${P_VALUE}) AS max_cpg,
                cpg
            FROM SNP_CPG_ARRAY
        )
        -- For 3 CpG per DMR, half of the SNPs are dropped
        SELECT * FROM HET_SNP WHERE nb_cpg >= ${CPG_PER_DMR}
        "

echo "Make a table with all the possible DMRs (before computing p-value)"
bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_snp_for_dmr \
    --replace=true \
    "
        -- Select SNPs where there are at least 3 significant CpGs in the same direction
        -- At this point, we're left with 2% of SNPs.
    WITH    
        HET_SNP AS (
            SELECT * FROM ${DATASET_ID}.${SAMPLE}_het_snp
        ),
        -- Any DMR needs at least CPG_PER_DMR in the same direction
        SNP_FOR_DMR AS (
            SELECT
                snp_id AS snp_id_dmr, 
                chr AS chr_dmr,
                min_cpg,
                max_cpg,
                nb_cpg,
                nb_sig_cpg,
                pos_sig_cpg,
                neg_sig_cpg,
                cpg
            FROM HET_SNP
            -- To have boundaries in the DMR made of significant CpGs, we need at least 2 of them.
            WHERE nb_sig_cpg >= 2 
        ),
        -- Import the list of CpGs with their respective snp_id, read_id, and allele
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
        WELL_COVERED_CPG AS (
            SELECT 
                chr AS chr_well_cov, 
                pos AS pos_well_cov
            FROM ${DATASET_ID}.${SAMPLE}_cpg_asm
        ),
        FILTERED_CPG AS (
            SELECT * FROM ALL_CPG
            INNER JOIN WELL_COVERED_CPG
            ON chr_cpg = chr_well_cov 
               AND pos_cpg = pos_well_cov 
        ),
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
        -- INTERSECT WITH THE LIST OF CPGS WHERE WE HAVE THE SNP AND READ_ID, AND ALLELE
        CPG_DMR AS (
            SELECT * FROM SNP_FOR_DMR
            INNER JOIN FILTERED_CPG_CLEAN
            ON snp_id_cpg = snp_id_dmr AND chr_cpg = chr_dmr
        ),
        QUALIFYING_CPG AS (
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
                pos_cpg >= min_cpg 
                AND pos_cpg <= max_cpg
        ),
        -- Compute the methylation % per read_id, per snp_id
        METHYL_PER_READ AS (
            SELECT 
                snp_id,
                chr_cpg,
                ANY_VALUE(min_cpg) AS dmr_inf,
                ANY_VALUE(max_cpg) AS dmr_sup,
                read_id,
                allele,
                ROUND(SAFE_DIVIDE(SUM(meth),SUM(cov)),5) AS methyl,
                ANY_VALUE(nb_cpg) AS nb_cpg,
                ANY_VALUE(nb_sig_cpg) AS nb_sig_cpg,
                ANY_VALUE(pos_sig_cpg) AS pos_sig_cpg,
                ANY_VALUE(neg_sig_cpg) AS neg_sig_cpg,
                ANY_VALUE(cpg) AS cpg
            FROM QUALIFYING_CPG
            GROUP BY snp_id, read_id, allele, chr_cpg
        ),
        SNP_METHYL_ARRAY_REF AS (
            SELECT
                snp_id,
                ANY_VALUE(chr_cpg) AS chr,
                ANY_VALUE(dmr_inf) AS dmr_inf,
                ANY_VALUE(dmr_sup) AS dmr_sup,
                ARRAY_AGG(STRUCT(methyl)) AS ref,
                ANY_VALUE(nb_cpg) AS nb_cpg,
                ANY_VALUE(nb_sig_cpg) AS nb_sig_cpg,
                ANY_VALUE(pos_sig_cpg) AS pos_sig_cpg,
                ANY_VALUE(neg_sig_cpg) AS neg_sig_cpg,
                ANY_VALUE(cpg) AS cpg
            FROM METHYL_PER_READ
            WHERE allele = 'REF'
            GROUP BY snp_id
        ),
        SNP_METHYL_ARRAY_ALT AS (
            SELECT
                snp_id AS snp_id_alt,
                ARRAY_AGG(STRUCT(methyl)) AS alt
            FROM METHYL_PER_READ
            WHERE allele = 'ALT'
            GROUP BY snp_id
        ),
        SNP_METHYL_JOIN AS (
            SELECT * FROM SNP_METHYL_ARRAY_REF
            INNER JOIN SNP_METHYL_ARRAY_ALT
            ON snp_id = snp_id_alt
        ),
        SNP_METHYL AS (
            SELECT 
                snp_id, 
                chr,
                dmr_inf,
                dmr_sup,
                ARRAY_LENGTH(ref) AS ref_reads, 
                ARRAY_LENGTH(alt) AS alt_reads,
                 ROUND(((SELECT AVG(methyl) FROM UNNEST(alt)) - (SELECT AVG(methyl) FROM UNNEST(ref))),3) AS effect,
                ref, 
                alt,
                nb_cpg,
                nb_sig_cpg,
                pos_sig_cpg,
                neg_sig_cpg,
                cpg
            FROM SNP_METHYL_JOIN
        )
        SELECT * FROM SNP_METHYL 
    "

# Export file to JSON format in the bucket
# (nested arrays are not supported in)
bq extract \
    --destination_format NEWLINE_DELIMITED_JSON \
    ${DATASET_ID}.${SAMPLE}_snp_for_dmr \
    gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_snp_for_dmr.json