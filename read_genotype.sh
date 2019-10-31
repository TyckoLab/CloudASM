#!/bin/bash

# CIGARs strings always start with a "M" (match). They are like 234M3D32M
# Then they may have an insertion or deletion and then a match

# Note: if the SNP is the latest nucleotide in the sequence, SUBSTR returns an empty value
# which we need to eliminate (only occurs in 0.02% of cases)

bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_vcf_reads_genotype \
    --replace=true \
    "
    WITH 
        -- transform the vcf_reads table by extracting the CIGAR string
        -- select only rows where the SNP is actually in the read. The other reads will be tagged later with REF or ALT based on their respective R1/R2
        VCF_ALL_READS AS (
            SELECT 
                snp_id,
                ref,
                alt,
                chr,
                pos,
                CT_strand,
                GA_strand,
                read_start,
                read_end,
                read_id,
                cigar,
                -- we split the CIGAR in 2 arrays, one with letters, one with their corresponding numbers.
                REGEXP_EXTRACT_ALL(cigar, r'[M|D|I]') AS a_cig_letters, 
                REGEXP_EXTRACT_ALL(cigar, r'[0-9]+') AS a_cig_num, 
                pos - read_start + 1 AS snp_pos_in_read,
                seq_CT_strand,
                seq_GA_strand,
                seq,
                score_before_recal
            FROM ${DATASET_ID}.${SAMPLE}_vcf_reads
            
        ),
        VCF_IN_READS AS (
            SELECT * FROM VCF_ALL_READS
            WHERE pos >= read_start AND pos <= read_end 
        ),
        CIG_NUM AS (
            SELECT 
                -- create a column for the number of letters in the CIGAR string (1, 3, 5, or 7)
                ARRAY_LENGTH(a_cig_num) AS length_cig,
                *
            FROM VCF_IN_READS
        ),
        CIG_NUM_REFINED AS (
            SELECT
                SAFE_CAST(a_cig_num[offset(0)] AS INT64) AS first_m,
                IF(length_cig >= 3,
                    SAFE_CAST(a_cig_num[offset(0)] AS INT64) 
                    + SAFE_CAST(a_cig_num[offset(2)] AS INT64)
                    + IF( a_cig_letters[offset(1)] = 'I', -1, 1)*SAFE_CAST(a_cig_num[offset(1)] AS INT64),
                    0)
                    AS second_m,
                IF(length_cig >= 5,
                            SAFE_CAST(a_cig_num[offset(0)] AS INT64)
                            + SAFE_CAST(a_cig_num[offset(2)] AS INT64)
                            + SAFE_CAST(a_cig_num[offset(4)] AS INT64)
                            + IF(a_cig_letters[offset(1)] = 'I', 1, -1) * SAFE_CAST(a_cig_num[offset(1)] AS INT64)
                            + IF(a_cig_letters[offset(3)] = 'I', 1, -1) * SAFE_CAST(a_cig_num[offset(3)] AS INT64),
                            0)
                    AS third_m,
            *
        FROM CIG_NUM
        ),
        -- first table: the SNP is located before the first deletion or insertion
        FIRST_MATCH AS (
            SELECT
                SUBSTR(seq, snp_pos_in_read, 1) as snp_in_read,
                -- We add 5 to take into account the tag OQ:Z:
                SUBSTR(score_before_recal, 5 + snp_pos_in_read, 1) AS snp_score,
                *
            FROM 
                CIG_NUM_REFINED
            WHERE
                first_m >= snp_pos_in_read
        ),
        -- second table: the SNP is located after the first deletion or insertion
        SECOND_MATCH AS (
            SELECT 
            SUBSTR(seq, snp_pos_in_read + IF( a_cig_letters[offset(1)] = 'I', 1, -1)*SAFE_CAST(a_cig_num[offset(1)] AS INT64), 1) as snp_in_read,
            SUBSTR(score_before_recal, 5 + snp_pos_in_read + IF( a_cig_letters[offset(1)] = 'I', 1, -1)*SAFE_CAST(a_cig_num[offset(1)] AS INT64), 1) as snp_score,
            *
            FROM 
                CIG_NUM_REFINED
            WHERE
                first_m < snp_pos_in_read AND second_m >= snp_pos_in_read
        ),
        -- third table: the SNP is located after the second deletion or insertion
        THIRD_MATCH AS (
            SELECT
                SUBSTR(seq, snp_pos_in_read + IF( a_cig_letters[offset(1)] = 'I', 1, -1)*SAFE_CAST(a_cig_num[offset(1)] AS INT64) + IF( a_cig_letters[offset(3)] = 'I', 1, -1)*SAFE_CAST(a_cig_num[offset(3)] AS INT64), 1) as snp_in_read,
                SUBSTR(score_before_recal, 5+ snp_pos_in_read + IF( a_cig_letters[offset(1)] = 'I', 1, -1)*SAFE_CAST(a_cig_num[offset(1)] AS INT64) + IF( a_cig_letters[offset(3)] = 'I', 1, -1)*SAFE_CAST(a_cig_num[offset(3)] AS INT64), 1) as snp_score,
                *
            FROM 
                CIG_NUM_REFINED
            WHERE
                first_m < snp_pos_in_read AND second_m < snp_pos_in_read AND third_m >= snp_pos_in_read
        ),
    FOR_GENOTYPING AS (
    SELECT 
        snp_id,
        snp_score,
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
    WHERE BYTE_LENGTH(snp_in_read) = 1 
    UNION ALL 
    SELECT 
        snp_id,
        snp_score,
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
    -- We demand that a letter was found
    WHERE BYTE_LENGTH(snp_in_read) = 1 
    UNION ALL 
    SELECT 
        snp_id,
        snp_score,
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
    -- We demand that a letter was found
    WHERE BYTE_LENGTH(snp_in_read) = 1 
    ),
    -- table when the SNP can be found on only the CT or GA strands. We require that the SNP found in the read is indeeed equal to REF or ALT. In some rare cases, it is not.
    ONLY_CT_OR_GA AS ( 
      SELECT
        snp_id,
        snp_score,
        chr,
        pos,
        read_id,
        if (snp_in_read = ref, 'REF', if (snp_in_read = alt, 'ALT', 'bad_snp')) AS allele
      FROM 
        FOR_GENOTYPING 
      WHERE 
       CT_strand = FALSE OR GA_STRAND = FALSE 
      ),
    BOTH_STRANDS AS (
    -- table when the SNP can be found on both strands
        SELECT 
            snp_id,
            snp_score,
            chr,
            pos,
            read_id,
            if (snp_in_read = ref 
                OR (ref = 'C' AND seq_CT_strand = TRUE AND snp_in_read = 'T') 
                OR (ref = 'G' AND seq_GA_strand = TRUE AND snp_in_read = 'A'), 
                    'REF', 
                        if (snp_in_read = alt
                        OR (alt = 'C' AND seq_CT_strand = TRUE AND snp_in_read = 'T')
                        OR (alt = 'G' AND seq_GA_strand = TRUE AND snp_in_read = 'A')
                            , 'ALT','bad_snp')) AS allele
        FROM FOR_GENOTYPING 
        WHERE 
            CT_strand = TRUE AND GA_STRAND = TRUE
    ),
    -- remove the SNPs that could be not be identified for sure.
    ALL_SNP AS (
        SELECT * FROM BOTH_STRANDS WHERE allele != 'bad_snp'
        UNION ALL
        SELECT * FROM ONLY_CT_OR_GA WHERE allele != 'bad_snp'
    ),
    -- This filters out the reads where the nucleotide with the SNP has a score below $SNP_SCORE defined in the master script
    SNPS_WITH_GOOD_QUALITY AS (
        SELECT DISTINCT * FROM ALL_SNP
        WHERE SAFE_CAST(TO_CODE_POINTS(snp_score)[offset(0)] AS INT64) >= ${SNP_SCORE}
    ),
    SNP_AND_READ_GENOTYPE AS (
        SELECT snp_id AS snp_id_tmp, read_id AS read_id_tmp, allele
        FROM SNPS_WITH_GOOD_QUALITY
    ),
    SNP_AND_ALL_READS AS (
        SELECT snp_id, chr, pos, read_id, allele FROM SNP_AND_READ_GENOTYPE
        INNER JOIN VCF_ALL_READS
        ON snp_id_tmp = snp_id AND read_id_tmp = read_id
    ),
    -- Need to remove the (snp_id, read_id) combination where the snp_id has been found to be both ref and alt (that can happen with paired reads both overlapping the snp position)
    GROUP_BY_READID AS (
        SELECT snp_id AS snp_id_rm, COUNT(DISTINCT allele) AS ambiguity, ANY_VALUE(read_id) AS read_id_rm
        FROM SNP_AND_ALL_READS
        GROUP BY read_id, snp_id
    ),
    SNP_READ_TO_EXCLUDE AS (
        SELECT snp_id_rm AS snp_id, read_id_rm AS read_id 
        FROM GROUP_BY_READID
        WHERE ambiguity > 1
    ),
    SNP_READ_TO_KEEP AS (
    SELECT snp_id, read_id FROM SNP_AND_ALL_READS
    EXCEPT DISTINCT
    SELECT snp_id, read_id FROM SNP_READ_TO_EXCLUDE
    ),
    REFORMAT AS (
    SELECT snp_id AS snp_id_keep, read_id AS read_id_keep 
    FROM SNP_READ_TO_KEEP
    ),
    FINAL_LIST_WITH_DUPLICATES AS (
    SELECT snp_id, chr, pos, read_id, allele FROM REFORMAT 
    INNER JOIN SNP_AND_ALL_READS
    ON snp_id = snp_id_keep AND read_id = read_id_keep
    )
    SELECT DISTINCT * FROM FINAL_LIST_WITH_DUPLICATES
  "



