#!/bin/bash

# Make sure it is consistent with what is created in the Docker image.
JAVA="/genomics-packages/jre1.7.0_25/bin/"
BIS_SNP="/genomics-packages/BisSNP-0.82.2"

# Temp folder for the packages
TMP_DIR="/mnt/data/tmp/"
mkdir -p ${TMP_DIR}

############## Variant call

# ‘-mmq 30’ specifies that reads with minimum mapping quality score more 
# than 30 would be used for genotyping.

# ‘-mbq 0’ specifies that the bases with minimum base quality score more than 0  
# are used for genotyping and methylation calling

# -nt is the number of threads.

# Variant call (mmq is the quality score)
$JAVA/java -Xmx48g \
    -Djava.io.tmpdir=${TMP_DIR} \
    -jar ${BIS_SNP}/BisSNP-0.82.2.jar \
    -L ${CHR} \
    -T BisulfiteGenotyper \
    -R $(dirname "${REF_GENOME}")/human_g1k_v37.fasta \
    -D ${VCF} \
    -I $(dirname "${BAM_BAI}")/${SAMPLE}_chr${CHR}_recal.bam \
    -vfn1 $(dirname "${BAM_BAI}")/${SAMPLE}.raw_chr$CHR.vcf \
    -mmq 30 \
    -mbq 0 \
    -stand_call_conf 20 \
    -nt 10 \
    -out_modes EMIT_HET_SNPS_ONLY

# Sorting the VCF by coordinate..."
perl ${BIS_SNP}/sortByRefAndCor.pl \
    --k 1 \
    --c 2 \
    --tmp ${TMP_DIR} \
    $(dirname "${BAM_BAI}")/${SAMPLE}.raw_chr$CHR.vcf \
    $(dirname "${REF_GENOME}")/human_g1k_v37.fasta.fai \
    > $(dirname "${BAM_BAI}")/${SAMPLE}.raw.sort_chr$CHR.vcf

# maxCOV: beyond this depth value, we ignore because we think it's a wrong region.

# Remove false positives by removing the calls on super high-coverage regions"
$JAVA/java  -Xmx48g -Djava.io.tmpdir=${TMP_DIR} \
    -jar ${BIS_SNP}/BisSNP-0.82.2.jar \
    -L $CHR \
    -R $(dirname "${REF_GENOME}")/human_g1k_v37.fasta \
    -T VCFpostprocess \
    -oldVcf $(dirname "${BAM_BAI}")/${SAMPLE}.raw.sort_chr$CHR.vcf \
    -newVcf $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_allvariants.vcf \
    -snpVcf $(dirname "${BAM_BAI}")/${SAMPLE}.raw.sort_chr$CHR.vcf \
    -o $(dirname "${OUTPUT_DIR}")/${SAMPLE}.snp.summary_chr$CHR.txt \
    -maxCov 200 \
    -minSNPinWind 2

# Import the VCF file in Big Query
bq --location=US load \
               --replace \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 117 \
               ${DATASET_ID}.${SAMPLE}_chr${CHR}_vcf \
               gs://$BUCKET/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}.vcf \
               chr:STRING,pos:STRING,snp_id:STRING,ref:STRING,alt:STRING,qual:FLOAT,filter:STRING,info:STRING,format:STRING,data:STRING

# Create a 500bp window around each SNP
bq query \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_${CHR}_500bp \
    --replace=true \
    "SELECT chr, 
       CAST(pos as INT64) - 250 as inf,
       CAST(pos as INT64) + 250 as sup,
       snp_id,
       ref,
       alt,
       CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) as ref_n,
       CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) as alt_n,
       CONCAT(chr,':', pos,'-',pos) as coord
    FROM `hackensack-tyco.wgbs_asm.varianttest`
    WHERE 
        snp_id LIKE '%rs%'
        AND CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) >= 10"


# Export 500bp window around each SNP to the bucket
bq extract \
    --field_delimiter "\t" \
    --print_header=false \
    ${DATASET_ID}.${SAMPLE}_${CHR}_500bp \
    gs://$OUTPUT_B/$SAMPLE/variants_per_chr/${SAMPLE}_${CHR}_500bp.bed
