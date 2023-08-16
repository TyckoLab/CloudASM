#!/bin/bash

# Sort BAM file (required before using samtools view)
samtools sort \
    $BAM \
    $(dirname ${BAM})/$(basename ${BAM%.bam})'_sorted'

# Create index file (required to split by chromosome)
samtools index $(dirname ${BAM})/$(basename ${BAM%.bam})'_sorted.bam'

########## Split the BAM file per chromosome

# Need to create these variables to avoid "ambiguous redirects" in the for loop
SORTED_BAM=$(dirname ${BAM})/$(basename ${BAM%.bam})'_sorted.bam'
FOLDER=$(dirname ${OUTPUT_DIR})
FASTQ=$(basename ${BAM%.bam})

# Check for the presence of 'chr' in SORTED_BAM
if samtools view ${SORTED_BAM} | grep -q "chr"; then
    # 'chr' is present in the BAM file
    CHR_PREFIX="chr"
else
    # 'chr' is not present in the BAM file
    CHR_PREFIX=""
fi

# Loop on chromosomes
for CHR in `seq 1 22` X Y ; do
    samtools view -bh ${SORTED_BAM} ${CHR_PREFIX}${CHR} > ${FOLDER}/${FASTQ}_chr${CHR}.bam
done

