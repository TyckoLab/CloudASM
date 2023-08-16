#!/bin/bash

# Note: Picardtools packages require that we specify the whole path of the package
# Make sure it is consistent with what is created in the Docker image.
PICARD="/genomics-packages/picard-tools-1.97"
JAVA="/genomics-packages/jre1.7.0_25/bin/"
BIS_SNP="/genomics-packages/BisSNP-0.82.2"

# Temp folder for the packages
TMP_DIR="/mnt/data/tmp/"
mkdir -p ${TMP_DIR}

# Sort by coordinate (required by Bis-SNP)
java \
      -Djava.io.tmpdir=${TMP_DIR} \
      -Xmx52g \
      -jar $PICARD/SortSam.jar \
      I=${BAM} \
      O=$(dirname "${BAM}")/${SAMPLE}_chr${CHR}_sorted.bam \
      SORT_ORDER=coordinate \
      CREATE_INDEX=true \
      MAX_RECORDS_IN_RAM=9000000

echo "Create recal file. This requires a specific version of Java."
java \
    -Djava.io.tmpdir=${TMP_DIR} \
    -Xmx52g \
    -jar ${BIS_SNP}/BisSNP-0.82.2.jar \
    -L chr$CHR \
    -R $(dirname "${REF_GENOME}")/*.fasta \
    -I $(dirname "${BAM}")/${SAMPLE}_chr${CHR}_sorted.bam \
    -T BisulfiteCountCovariates \
    -knownSites ${ALL_VARIANTS} \
    -cov ReadGroupCovariate \
    -cov QualityScoreCovariate \
    -cov CycleCovariate \
    -recalFile $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_recal.csv \
    -nt 12

echo  "Generates a re-calibrated BAM and BAI file"
java \
    -Djava.io.tmpdir=${TMP_DIR} \
    -Xmx52g \
    -jar ${BIS_SNP}/BisSNP-0.82.2.jar \
    -L chr$CHR \
    -R $(dirname "${REF_GENOME}")/*.fasta \
    -I $(dirname "${BAM}")/${SAMPLE}_chr${CHR}_sorted.bam \
    -o $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_recal.bam \
    -T BisulfiteTableRecalibration \
    -recalFile $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_recal.csv \
    --max_quality_score 40 
