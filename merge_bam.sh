#!/bin/bash

# Note: Picardtools packages require that we specify the whole path of the package
# Make sure it is consistent with what is created in the Docker image.
PICARD="/genomics-packages/picard-tools-1.97"

# Define temporary directory
TMP_DIR="/mnt/data/tmp/"
mkdir -p ${TMP_DIR}

echo "Merge per chromosome..."
NB_BAM_FILES=$(ls ${BAM_FILES} | wc -l)
if [ ${NB_BAM_FILES} == 1 ] ; then
    echo "Only one chard is relevant"
    cp ${BAM_FILES} $(dirname "${TMP_DIR}")/${SAMPLE}_chr${CHR}_merged.bam
else 
    echo "Several chards are relevant"
    samtools merge $(dirname "${TMP_DIR}")/${SAMPLE}_chr${CHR}_merged.bam $BAM_FILES
fi

echo "Sort by coordinate"
java -Xmx48g \
      -Djava.io.tmpdir=${TMP_DIR} \
      -jar $PICARD/SortSam.jar \
      I=$(dirname "${TMP_DIR}")/${SAMPLE}_chr${CHR}_merged.bam \
      O=$(dirname "${TMP_DIR}")/${SAMPLE}_chr${CHR}_sort.bam \
      SORT_ORDER=coordinate \
      MAX_RECORDS_IN_RAM=9000000

echo "Perform AddReplaceGroups (required by snp-calling software Bis-SNP)"
java -Xmx48g \
    -Djava.io.tmpdir=${TMP_DIR} \
    -jar $PICARD/AddOrReplaceReadGroups.jar \
    I=$(dirname "${TMP_DIR}")/${SAMPLE}_chr${CHR}_sort.bam \
    O=$(dirname "${TMP_DIR}")/${SAMPLE}_chr${CHR}_sort_ARG.bam \
    ID=${SAMPLE}_ID \
    LB=${SAMPLE}_LB \
    PL=illumina \
    PU=run \
    SM=$SAMPLE \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    SORT_ORDER=coordinate \
    MAX_RECORDS_IN_RAM=9000000

echo "Remove PCR duplicates"
java -Xmx48g \
    -Djava.io.tmpdir=${TMP_DIR} \
    -jar $PICARD/MarkDuplicates.jar \
    INPUT=$(dirname "${TMP_DIR}")/${SAMPLE}_chr${CHR}_sort_ARG.bam \
    OUTPUT=$(dirname "${TMP_DIR}")//${SAMPLE}_chr${CHR}_sort_ARG_RD.bam \
    METRICS_FILE=$(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_duplicates.txt \
    TMP_DIR=${TMP_DIR} \
    CREATE_INDEX=true \
    REMOVE_DUPLICATES=true \
    MAX_RECORDS_IN_RAM=9000000

echo "BAM files merged by chr before sorting by read ID"
ls -lh $(dirname "${TMP_DIR}")/

echo "Sort by read ID (option -n), required for methylation extractor"
samtools sort \
    -@ 14 \
    -n \
    $(dirname "${TMP_DIR}")/${SAMPLE}_chr${CHR}_sort_ARG_RD.bam \
    $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}
