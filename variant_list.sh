#!/bin/bash

# # Intersect the variant list with CpG positions
# intersectBed \
#         -a "${VCF_500BP}" \
#         -b "${CPG_POS}" \
#         > $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_500bp_wCpG.bed

# # Create a list of all the variants on which to call ASM
# # Columns are snp_id, chr, position, coverage, and coordinates with a 2000 window(chr:inf-sup)
# awk 'BEGIN { FS=OFS="\t" }
#         {print $4,$1,$8,$5,$6,$7,$9}' \
#         $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_500bp_wCpG.bed \
#          > $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_variants.txt


# Divide the list of SNPs into equal pieces of 200

split -l 200 \
        --numeric-suffixes --suffix-length=6 \
        --additional-suffix=.txt \
        $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_variants.txt \
        $(dirname "${OUTPUT_DIR}")/${SAMPLE}_chr${CHR}_variants_

