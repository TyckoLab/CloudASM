DOCKER_GENOMICS="gcr.io/hackensack-tyco/wgbs-asm"
DOCKER_R="gcr.io/hackensack-tyco/r-genomics"
DOCKER_PERL="gcr.io/hackensack-tyco/perl-samtools"

SCRIPTS="$HOME/GITHUB_REPOS/wgbs-asm/onuchic_scripts"

PROJECT_ID="hackensack-tyco"
REGION_ID="us-central1"
ZONE_ID="us-central1-b"
OUTPUT_B="em-encode-paper" # will be created by the script
LOG="gs://$OUTPUT_B/logging"
REF_DATA_B="wgbs-ref-files" # will be created by the script
GENOME="hg19" 

SAMPLE="gm12878mod"
CHR="1"

REF_GENOME="gs://em-encode-paper/gm12878mod/onuchic_ref_genome/human_g1k_v37_chr1.fasta"
REF_GENOME_INDEX="gs://em-encode-paper/gm12878mod/onuchic_ref_genome/human_g1k_v37_chr1.fasta.fai"


# Link non bisulfite-converted VCF for Hg19 for the gm12878 cell line
ftp://platgene_ro@ussd-ftp.illumina.com/2016-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz

# Sort BAM file 
dsub \
  --provider google-v2 \
  --ssh \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --logging $LOG \
  --disk-size 200 \
  --machine-type n1-standard-16 \
  --image=$DOCKER_GENOMICS \
  --input BAM="gs://$OUTPUT_B/gm12878/bam_per_chr/gm12878_chr1.bam" \
  --output BAM_SORTED="gs://$OUTPUT_B/gm12878mod/onuchic_bam_sorted/gm12878mod_chr1_sorted.bam" \
  --output BAM_INDEX="gs://$OUTPUT_B/gm12878mod/onuchic_bam_sorted/gm12878mod_chr1_sorted.bam.bai" \
  --script $SCRIPTS/bam_sort.sh \
  --wait



#----------------------------------------------------------------------------------------------------------------------------------

# Generate a list of CpG sites
# 10 min total. --> 1 CPU-hours
# 15h50
dsub \
  --provider google-v2 \
  --ssh \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --logging $LOG \
  --image=$DOCKER_PERL \
  --machine-type n1-standard-8 \
  --boot-disk-size 100 \
  --disk-size 100 \
  --input VCF="gs://em-encode-paper/gm12878mod/onuchic_vcf/NA12878_chr1_final.vcf" \
  --input REF_GENOME_INDEX="${REF_GENOME_INDEX}" \
  --input REF_GENOME="${REF_GENOME}" \
  --output-recursive OUTPUT_DIR="gs://$OUTPUT_B/$SAMPLE/onuchic_cpg" \
  --script $SCRIPTS/getCpgPositions.pl \
  --wait

# Download the HetSNp.txt file on a computer. Replace "chr" with nothing
sed 's|chr||g' hetSnps.txt > hetSnps_nochr.txt 
# Then, upload to bucket

# Methylation count per CpG
# Took one hour --> 8 CPU-hours
# 15h15
# The output is actually the stdout file -- it does not save any file. Upload the file as gs://$OUTPUT_B/$SAMPLE/onuchic_methyl/$SAMPLE_chr$CHR_methyl.txt
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
  --input BAM="gs://em-encode-paper/gm12878mod/onuchic_bam_sorted/gm12878mod_chr1_rs2319164.bam" \
  --input BAM_INDEX="gs://em-encode-paper/gm12878mod/onuchic_bam_sorted/gm12878mod_chr1_rs2319164.bam.bai" \
  --input REF_GENOME_INDEX="${REF_GENOME_INDEX}" \
  --input REF_GENOME="${REF_GENOME}" \
  --input CPGPOS="gs://$OUTPUT_B/$SAMPLE/onuchic_cpg/homoCpGs.txt" \
  --input SNPS="gs://$OUTPUT_B/$SAMPLE/onuchic_cpg/hetSnps_nochr.txt" \
  --script $SCRIPTS/getAllelicMethCounts.pl \
  --wait


