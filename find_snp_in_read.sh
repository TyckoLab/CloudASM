#!/bin/bash

# CIGARs strings always start with a "M" (match). They are like 234M3D32M
# Then they may have an insertion or deletion and then a match

bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_reads_tmp_${CHR} \
    --replace=true \
    "
    WITH 
        VCF_reads AS (
            SELECT 
                snp_id,
                ref,
                alt,
                chr,
                pos,
                CT_strand,
                GA_strand,
                read_id,
                cigar,
                read_start,
                -- we split the CIGAR in 2 arrays, one with letters, one with their corresponding numbers.
                REGEXP_EXTRACT_ALL(cigar, r'[M|D|I]') AS a_cig_letters, 
                REGEXP_EXTRACT_ALL(cigar, r'[0-9]+') AS a_cig_num, 
                pos - read_start + 1 AS snp_pos_in_read,
                seq_CT_strand,
                seq_GA_strand,
                seq
            FROM `wgbs_asm.gm12878_vcf_reads_tmp`
        ),
        CIG_NUM AS (
            SELECT 
                -- create a column for the number of letters in the CIGAR string (1, 3, 5, or 7)
                ARRAY_LENGTH(a_cig_num) AS length_cig,
                *
            FROM VCF_READS
        ),
        CIG_NUM_REFINED AS (
            SELECT
                SAFE_CAST(a_cig_num[offset(0)] AS INT64) AS first_m,
                IF(length_cig >= 3,
                    SAFE_CAST(a_cig_num[offset(0)] AS INT64) 
                    + SAFE_CAST(a_cig_num[offset(2)] AS INT64)
                    + IF( a_cig_letters[offset(1)] = 'I', 1, -1)*SAFE_CAST(a_cig_num[offset(1)] AS INT64),
                    0)
                    AS second_m,
                IF(length_cig >= 5,
                            SAFE_CAST(a_cig_num[offset(0)] AS INT64)
                            + SAFE_CAST(a_cig_num[offset(2)] AS INT64)
                            + SAFE_CAST(a_cig_num[offset(4)] AS INT64)
                            + IF(a_cig_letters[offset(1)] = 'I', 1, - 1) * SAFE_CAST(a_cig_num[offset(1)] AS INT64)
                            + IF(a_cig_letters[offset(3)] = 'I', 1, -1) * SAFE_CAST(a_cig_num[offset(3)] AS INT64),
                            0)
                    AS third_m,
            *
        FROM CIG_NUM
        ),
        FIRST_MATCH AS (
            SELECT
                SUBSTR(seq, snp_pos_in_read, 1) as snp_in_read,
                *
            FROM 
                CIG_NUM_REFINED
            WHERE
                first_m >= snp_pos_in_read
        ),
        SECOND_MATCH AS (
            SELECT 
            SUBSTR(seq, snp_pos_in_read + IF( a_cig_letters[offset(1)] = 'I', 1, -1)*SAFE_CAST(a_cig_num[offset(1)] AS INT64), 1) as snp_in_read,
            *
            FROM 
                CIG_NUM_REFINED
            WHERE
                first_m < snp_pos_in_read AND second_m >= snp_pos_in_read
        ),
        THIRD_MATCH AS (
            SELECT
                SUBSTR(seq, snp_pos_in_read + IF( a_cig_letters[offset(1)] = 'I', 1, -1)*SAFE_CAST(a_cig_num[offset(1)] AS INT64) + IF( a_cig_letters[offset(3)] = 'I', 1, -1)*SAFE_CAST(a_cig_num[offset(3)] AS INT64), 1) as snp_in_read,
                *
            FROM 
                CIG_NUM_REFINED
            WHERE
                first_m < snp_pos_in_read AND second_m < snp_pos_in_read AND third_m >= snp_pos_in_read
        ),
    FOR_GENOTYPING AS (
    SELECT 
        snp_id,
        ref,
        alt,
        chr,
        pos,
        CT_strand,
        GA_strand,
        snp_in_read,
        seq_CT_strand,
        seq_GA_strand,
        read_id
    FROM THIRD_MATCH
    UNION ALL 
    SELECT 
        snp_id,
        ref,
        alt,
        chr,
        pos,
        CT_strand,
        GA_strand,
        snp_in_read,
        seq_CT_strand,
        seq_GA_strand,
        read_id
    FROM SECOND_MATCH
    UNION ALL 
    SELECT 
        snp_id,
        ref,
        alt,
        chr,
        pos,
        CT_strand,
        GA_strand,
        snp_in_read,
        seq_CT_strand,
        seq_GA_strand,
        read_id
    FROM FIRST_MATCH
    ),
    ONLY_CT_OR_GA AS ( 
      SELECT
        snp_id,
        chr,
        read_id,
        if (snp_in_read = ref, 'REF', 'ALT') AS allele
      FROM FOR_GENOTYPING WHERE 
        -- if the SNP can be found on only CT or GA, the SNP found must be equal to the ref or alt reported by bis-snp. In some rare cases there is a problem, which we eliminate.
       ((CT_strand = FALSE OR GA_STRAND = FALSE) AND (snp_in_read = ref OR snp_in_read = alt)) 
    )
    SELECT * FROM ONLY_CT_OR_GA
  "