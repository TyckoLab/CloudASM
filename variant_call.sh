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

echo "Variant call (mmq is the quality score)"
$JAVA/java -Xmx50g \
    -Djava.io.tmpdir=${TMP_DIR} \
    -jar ${BIS_SNP}/BisSNP-0.82.2.jar \
    -L chr${CHR} \
    -T BisulfiteGenotyper \
    -R $(dirname "${REF_GENOME}")/*.fasta \
    -D ${ALL_VARIANTS} \
    -I $(dirname "${BAM_BAI}")/${SAMPLE}_chr${CHR}_recal.bam \
    -vfn1 $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_raw.vcf \
    -mmq 30 \
    -mbq 0 \
    -stand_call_conf 20 \
    -nt 12 \
    -out_modes EMIT_HET_SNPS_ONLY

echo "Sorting the VCF by coordinate..."
perl ${BIS_SNP}/sortByRefAndCor.pl \
    --k 1 \
    --c 2 \
    --tmp ${TMP_DIR} \
    $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_raw.vcf \
    $(dirname "${REF_GENOME}")/*.fasta.fai \
    > $(dirname "${BAM_BAI}")/${SAMPLE}.raw.sort_chr$CHR.vcf

# maxCOV: beyond this depth value, we ignore because we think it's a wrong region.

echo "folder content after sorting the VCF"
ls -lh $(dirname "${OUTPUT_DIR}")

echo "BAM_BAI"
ls -lh $(dirname "${BAM_BAI}")


echo "Remove false positives by removing the calls on super high-coverage regions"
# This removes about 7% of the SNPs identified in "raw"
# Note: the option -nt makes the following command fail.
$JAVA/java -Xmx50g \
    -Djava.io.tmpdir=${TMP_DIR} \
    -jar ${BIS_SNP}/BisSNP-0.82.2.jar \
    -L chr$CHR \
    -R $(dirname "${REF_GENOME}")/*.fasta \
    -T VCFpostprocess \
    -oldVcf $(dirname "${BAM_BAI}")/${SAMPLE}.raw.sort_chr$CHR.vcf \
    -newVcf $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_filtered.vcf \
    -snpVcf $(dirname "${BAM_BAI}")/${SAMPLE}.raw.sort_chr$CHR.vcf \
    -o $(dirname "${OUTPUT_DIR}")/${SAMPLE}.snp.summary_chr$CHR.txt \
    -maxCov 200 \
    -minSNPinWind 2
    

