#!/bin/bash

# Intersect the variant list with CpG positions
intersectBed \
        -a "${VCF_500BP}" \
        -b "${CPG_POS}" \
        > $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_500bp_wCpG.bed

# Create a list of all the variants on which to call ASM
# Columns are chr, snp_id, and coordinates (chr:inf-sup)
awk 'BEGIN { FS=OFS="\t" }
        {print $1,$4,$9}' \
        $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_500bp_wCpG.bed \
         > $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_variants.txt

# Slice that list in 10 shards
echo -e "1,50\n2,50\n3,50\n4,50\n5,40\n6,40\n7,40\n8,40\n9,30\n10,30\n11,30\n12,30\n13,20\n14,20\n15,20\n16,20\n17,15\n18,15\n19,15\n20,15\n21,10\n22,10\nX,5\nY,5" > shard_counts.txt

# Determine the number of shards as a function of the chromosome
SHARD_COUNT=$(awk -v CHR="${CHR}" -F "," '{if($1==CHR) print $2}' shard_counts.txt)

# Split the variant list in shards
for ((SHARD = 1; SHARD < ${SHARD_COUNT}; SHARD++)); do
  awk "(NR%${SHARD_COUNT})==${SHARD}"'{ print $0 }' \
    $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_variants.txt \
    > $(dirname "${OUTPUT_DIR}")/split-${SHARD}-of-${SHARD_COUNT}_${SAMPLE}_chr${CHR}.txt
done

# Special situation not taken into account by the loop
awk "(NR%${SHARD_COUNT})==0"'{ print $0 }' \
    $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_variants.txt \
    > $(dirname "${OUTPUT_DIR}")/split-${SHARD_COUNT}-of-${SHARD_COUNT}_${SAMPLE}_chr${CHR}.txt

