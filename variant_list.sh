#!/bin/bash

# Keep the variants with CpGs in their 500bp viscinity
intersectBed \
    -a $(dirname $BED_500)/${SAMPLE}_chr${CHR}_vcf_500bp.bed \
    -b ${CPG_POS} \
    -wa \
    -wb \
    | awk 'BEGIN{FS=OFS="\t"} {print $11, $1,$3-250,$4,$5,$6,$7,$8,$9}' \
    |  uniq > $(dirname $BED_500)/${SAMPLE}_chr${CHR}_500bp.vcftemp

# Creates a specific file for ASM calling. (chromosome, SNP id, coordinates)
awk 'BEGIN{FS=OFS="\t"} {print $1,$3,$8}' \
    ${SAMPLE}_chr${CHR}.vcf \
    | sed 1d \
    > snp_pos_${SAMPLE}_chr${CHR}

# ${SAMPLE}_chr${CHR}.vcf is the final FILE.

# Slices the SNP pos file in 54 pieces.
awk -v sample="${SAMPLE}_chr${CHR}" '{
    a[NR] = $0
}
END {
    for (i = 1; i in a; ++i) {
        x = (i * 54 - 1) / NR + 1
        sub(/\..*$/, "", x)
        print a[i] >  x"_snp_pos_"sample
    }
}' snp_pos_${SAMPLE}_chr${CHR}

# Creates a list of the files
seq 1 `ls *_snp_pos_${SAMPLE}_chr${CHR} | wc -l` \
    > ${SAMPLE}_chr${CHR}_list

##


## 39,596
## 239,754

gm12878_chr22_20.vcf: 24,951
${SAMPLE}_chr${CHR}_vcf_500bp.bed : 24,951
${SAMPLE}_chr${CHR}_500bp.vcftemp (AFTER INTERSECTBED) : 213,709 rows?? (create a unique combination of window + SNP)