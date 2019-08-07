#!/bin/bash

# Note: Picardtools packages require that we specify the whole path of the package
# Make sure it is consistent with what is created in the Docker image.
PICARD="/genomics-packages/picard-tools-1.97"
JAVA="/genomics-packages/jre1.7.0_25/bin/"

# Temp folder for the packages
TMP_DIR="/mnt/data/tmp/"
mkdir -p ${TMP_DIR}

# Sort by coordinate (required by Bis-SNP)
java -Xmx48g -jar $PICARD/SortSam.jar \
      I=${BAM} \
      O=$(dirname ${BAM})/${SAMPLE}_chr${CHR}_sorted.bam \
      SORT_ORDER=coordinate \
      MAX_RECORDS_IN_RAM=9000000

# Create index file
samtools index $(dirname ${BAM})/${SAMPLE}_chr${CHR}_sorted.bam

# Create recal file
$JAVA/java \
    -Djava.io.tmpdir=${TMP_DIR} \
    -jar -Xmx48g $PATH_BISNP/BisSNP-0.82.2.jar \
    -L $CHR \
    -R $REF_GENOME/human_g1k_v37.fasta \
    -I $(dirname ${BAM})/${SAMPLE}_chr${CHR}_sorted.bam \
    -T BisulfiteCountCovariates \
    -knownSites $PATH_VCF \
    -cov ReadGroupCovariate \
    -cov QualityScoreCovariate \
    -cov CycleCovariate \
    -recalFile $(dirname ${OUTPUT_DIR})/${SAMPLE}_chr${CHR}_recal.csv \
    -nt 10

# # Generates a re-calibrated BAM and BAI file
# $JAVA/java -Djava.io.tmpdir=${TMP_DIR} \
#     -jar -Xmx48g $PATH_BISNP/BisSNP-0.82.2.jar \
#     -L $CHR \
#     -R $REF_GENOME/human_g1k_v37.fasta \
#     -I $WD/downloaded/${SAMPLE}_chr${CHR}_sorted.bam \
#     -o ${SAMPLE}_merged_dedup_RR_recal_chr${CHR}.bam \
#     -T BisulfiteTableRecalibration \
#     -recalFile ${SAMPLE}_merged_dedup_RR_chr${CHR}.recal.csv \
#     --max_quality_score 40

# echo
# echo $(date) "Test if the recalibration worked"
# echo
# ## Creates several PDF of graphs
# $JAVA/java -Xmx48g \
#     -jar $PATH_BISNP/BisulfiteAnalyzeCovariates-0.69.jar \
#     -recalFile ${SAMPLE}_merged_dedup_RR_chr${CHR}.recal.csv \
#     -outputDir $WD \
#     -ignoreQ 5 \
#     --max_quality_score 40


# # Print folder 
# echo
# echo "************** CHECK THE FILES AFTER THE RECAL STEP *****************"
# ls -ahl $WD
# echo "*****************************************************"


# #------------------------------------------------------------------------
# # Copy data over to the bucket and the $HOME folder
# #------------------------------------------------------------------------

# echo
# echo "Uploading the output of the recal algo (recal BAM, recal BAI, and PDF reports)"
# gsutil -o GSUtil:parallel_composite_upload_threshold=150M \
#     cp $WD/*recal*.bam gs://$BUCKET/$SAMPLE/"merged"

# # Copy the index file
# gsutil cp $WD/*recal*.bai gs://$BUCKET/$SAMPLE/"merged"

# # Copy the csv file
# gsutil cp $WD/*recal.csv gs://$BUCKET/$SAMPLE/"merged"


# # Copy the PDF reports
# gsutil -m cp $WD/*pdf  gs://$BUCKET/$SAMPLE/"recal_reports"/${CHR}

# # Copy the dat file
# gsutil cp $WD/*.dat gs://$BUCKET/$SAMPLE/"recal_reports"/${CHR}

# gsutil 

# #------------------------------------------------------------------------
# # Print message at the end
# #------------------------------------------------------------------------

# # Space left on the local disk
# echo
# echo "**************    DISK LIST AFTER JOB    *****************"
# df -h
# echo "**********************************************************"


# echo
# echo "End of script to create the recal BAM file"
# echo "Job finished on" $(date)

