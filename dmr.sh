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
        )
       -- SELECT snp_id AS snp_id_for_dmr, chr, nb_cpg, min_cpg, max_cpg FROM DMR_LIST WHERE pos_sig_cpg >=3 OR neg_sig_cpg >=3
        SELECT * FROM DMR_LIST
    "

# WITH TEST AS (
#   SELECT  *  FROM `hackensack-tyco.wgbs_asm.gm12878_cpg_asm` 
#   )
#  SELECT TO_JSON_STRING(TEST) FROM TEST



Wilcoxon from ASM script


        # Create an array with chr #, snp ID, first cpg, last cpg
        summary.snp = byread.snp[!duplicated(byread.snp$snp_id), c(2, 6, 4:5)]

        # Add two columns with the average methylation per haplotype 
        summary.snp$avg_meth_read_ref = avg.methyl.per.hapl[which(avg.methyl.per.hapl$haplotype == "ref"), 1]
        summary.snp$avg_meth_read_alt = avg.methyl.per.hapl[which(avg.methyl.per.hapl$haplotype == "alt"), 1]
        
        # Calculate the difference in the average methylation per haplotype (ALT - REF)
        summary.snp$diff_alt_ref = summary.snp$avg_meth_read_alt - summary.snp$avg_meth_read_ref
        
        summary.snp$wilcox = wilcox.test(meth_per_read ~ haplotype, data = byread.snp)$p.value

Compute av. methylation (not av. methylation %) across all CpGs between the 2 significant one