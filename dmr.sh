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

# Query to select the SNPs with at least CPG_PER_DMR nearby.
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
                ARRAY_AGG(STRUCT(pos,
                                ROUND(alt_meth/alt_cov-ref_meth/ref_cov,3) AS alt_minus_ref,
                                fisher_pvalue)) 
                        AS cpg
            FROM ${DATASET_ID}.${SAMPLE}_cpg_asm
            GROUP BY snp_id
        ),
        -- Extract the number of CpG per SNP 
        HET_SNP AS (
            SELECT 
                snp_id, 
                chr, 
                ARRAY_LENGTH(cpg) AS nb_cpg, 
                (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg) WHERE fisher_pvalue < 0.05) AS nb_sig_cpg, 
                (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg) WHERE fisher_pvalue < 0.05 AND SIGN(alt_minus_ref) = 1) AS pos_sig_cpg,
                (SELECT COUNT(fisher_pvalue) FROM UNNEST(cpg) WHERE fisher_pvalue < 0.05 AND SIGN(alt_minus_ref) = -1) AS neg_sig_cpg, 
                (SELECT MIN(pos) FROM UNNEST(cpg) WHERE fisher_pvalue < 0.05) AS min_cpg,
                (SELECT MAX(pos) FROM UNNEST(cpg) WHERE fisher_pvalue < 0.05) AS max_cpg
            FROM SNP_CPG_ARRAY
        )
        -- For 3 CpG per DMR, half of the SNPs are dropped
        SELECT * FROM HET_SNP WHERE nb_cpg >= ${CPG_PER_DMR}
        "


bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_snp_for_dmr \
    --replace=true \
    "
        -- Select SNPs where there are at least 3 significant CpGs in the same direction
        -- At this point, we're left with 2% of SNPs.
        SNP_FOR_DMR AS (
            SELECT
                snp_id AS snp_id_dmr, 
                chr AS chr_dmr,
                min_cpg,
                max_cpg
            FROM SNP_DETAILS
            WHERE 
                pos_sig_cpg >= ${CPG_PER_DMR}
                OR neg_sig_cpg >= ${CPG_PER_DMR}
        ),
        -- Import the list of CpGs with their respective snp_id, read_id, and allele
        CPG_LIST AS (
            SELECT * FROM ${DATASET_ID}.${SAMPLE}_cpg_read_genotype
        ),
        -- INTERSECT WITH THE LIST OF CPGS WHERE WE HAVE THE SNP AND READ_ID, AND ALLELE
        SNP_AND_READID AS (
            SELECT * FROM SNP_FOR_DMR
            INNER JOIN CPG_LIST
            ON snp_id = snp_id_dmr AND chr = chr_dmr
        ),
        QUALIFYING_CPG AS (
        -- All CpGs that are located within the boundaries of the potential DMR
        -- This removes 1/3 of all CpGs.
            SELECT 
                chr, 
                pos, 
                meth, 
                cov, 
                snp_id, 
                allele, 
                read_id 
            FROM SNP_AND_READID
            WHERE 
                pos >= min_cpg 
                AND pos <= max_cpg
        ),
        -- Compute the methylation % per read_id, per snp_id
        METHYL_PER_READ AS (
            SELECT 
                snp_id,
                chr,
                read_id,
                allele,
                ROUND(SAFE_DIVIDE(SUM(meth),SUM(cov)),5) AS methyl
            FROM QUALIFYING_CPG
            GROUP BY snp_id, read_id, allele, chr
        ),
        SNP_METHYL_ARRAY_REF AS (
            SELECT
                snp_id,
                ANY_VALUE(chr) AS chr,
                ARRAY_AGG(STRUCT(methyl)) AS ref
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
                ARRAY_LENGTH(ref) AS ref_reads, 
                ARRAY_LENGTH(alt) AS alt_reads,
                 ROUND(((SELECT AVG(methyl) FROM UNNEST(alt)) - (SELECT AVG(methyl) FROM UNNEST(ref))),3) AS effect,
                ref, 
                alt
            FROM SNP_METHYL_JOIN
        )
        -- This removes about 15% of potential DMR
        --SELECT * FROM SNP_METHYL WHERE abs(effect) > ${DMR_EFFECT}
        SELECT * FROM METHYL_PER_READ 
    "

# Export file to JSON format in the bucket
# (nested arrays are not supported in)
bq extract \
    --destination_format NEWLINE_DELIMITED_JSON \
    ${DATASET_ID}.${SAMPLE}_snp_for_dmr \
    gs://$OUTPUT_B/$SAMPLE/asm/${SAMPLE}_snp_for_dmr.json