DOCKER_R="gcr.io/hackensack-tyco/r-genomics"
DOCKER_PERL="gcr.io/hackensack-tyco/r-genomics"


# Download the scripts
git clone https://github.com/BRL-BCM/allelic_epigenome.git

ONUCHIC=$HOME/"GITHUB_REPOS"/"allelic_epigenome"

# Generate a list of CpG sites
# 20 min total.
dsub \
  --provider google-v2 \
  --ssh \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --logging $LOG \
  --image=$DOCKER_PERL \
  --machine-type n1-standard-8 \
  --boot-disk-size 200 \
  --disk-size 200 \
  --input VCF="gs://$OUTPUT_B/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}_filtered.vcf" \
  --input REFGEN="${REF_GENOME}/human_g1k_v37.fasta" \
  --output-recursive OUTPUT_DIR="gs://$OUTPUT_B/$SAMPLE/onuchic_cpg" \
  --script $ONUCHIC/Perl/getCpgPositions.pl \
  --wait


dsub \
  --provider google-v2 \
  --ssh \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --logging $LOG \
  --image=$DOCKER_PERL \
  --machine-type n1-standard-8 \
  --boot-disk-size 200 \
  --disk-size 200 \
  --input BAM="gs://$OUTPUT_B/$SAMPLE/recal_bam_per_chr/${SAMPLE}_chr${CHR}_recal.bam" \
  --input REFGEN="${REF_GENOME}/human_g1k_v37.fasta" \
  --input CPGPOS="gs://$OUTPUT_B/$SAMPLE/onuchic_cpg/homoCpGs.txt" \
  --input SNPS="gs://$OUTPUT_B/$SAMPLE/onuchic_cpg/hetSnps.txt" \
  --output-recursive OUTPUT_DIR="gs://$OUTPUT_B/$SAMPLE/onuchic_methyl" \
  --script $ONUCHIC/Perl/getAllelicMethCounts.pl \
  --wait


# Locally on a machine.
VCF="/wgbs/gm12878mod_chr1_filtered.vcf"
REF_GENOME="/wgbs/ref_genome/human_g1k_v37.fasta"
perl getCpgPositions.pl $REF_GENOME $VCF "gm12878chr1"