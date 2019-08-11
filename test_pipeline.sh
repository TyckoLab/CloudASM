########################## Simple test for dsub ################################

# This should not lead to any error
dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --regions $REGION_ID \
  --image $DOCKER_IMAGE \
  --logging gs://$OUTPUT_B/logging/ \
  --output OUT=gs://$OUTPUT_B/output/out.txt \
  --command 'echo "Hello World" > "${OUT}"' \
  --wait

# Create a folder with write permissions for the logging
mkdir "/tmp"
export TMPDIR="/tmp"
chmod -R 755 $TMPDIR

# If the logging occurs online, add /var/folders to Docker paths.

dsub \
  --provider local \
  --logging gs://$OUTPUT_B/logging/ \
  --command 'echo "Hello World"' \
  --wait

# Delete the output
gsutil rm -r gs://$OUTPUT_B/output



########################## Make an extract on ENCODE to test the pipeline ################################

# Type this in extract.tsv

--input ZIPPED	--output EXTRACT
gs://encode-wgbs/A549/ENCFF327QCK.fastq.gz	gs://encode-wgbs/A549-extract/ENCFF327QCK.fastq.gz
gs://encode-wgbs/A549/ENCFF986UWM.fastq.gz	gs://encode-wgbs/A549-extract/ENCFF986UWM.fastq.gz
gs://encode-wgbs/gm12878/ENCFF798RSS.fastq.gz	gs://encode-wgbs/gm12878-extract/ENCFF798RSS.fastq.gz
gs://encode-wgbs/gm12878/ENCFF113KRQ.fastq.gz	gs://encode-wgbs/gm12878-extract/ENCFF113KRQ.fastq.gz

dsub \
  --provider google-v2 \
  --project $PROJECT_ID \
  --zones $ZONE_ID \
  --logging gs://$OUTPUT_B/logging/ \
  --disk-size 800 \
  --preemptible \
  --image $DOCKER_IMAGE \
  --command 'gunzip ${ZIPPED} && \
             head -12000000 ${ZIPPED%.gz} > ${EXTRACT%.gz} && \
             gzip ${EXTRACT%.gz}' \
  --tasks extract.tsv \
  --wait







  ########################## TRIM GALORE ON ONE EXAMPLE

  # dsub \
#   --provider google-v2 \
#   --project $PROJECT_ID \
#   --zones $ZONE_ID \
#   --image $DOCKER_IMAGE \
#   --machine-type n1-standard-2 \
#   --logging gs://$OUTPUT_B/logging/ \
#   --input R1="gs://$OUTPUT_B/A549-extract/split_fastq/A549-extract_L01.R1.001.fastq" \
#   --input R2="gs://$OUTPUT_B/A549-extract/split_fastq/A549-extract_L01.R2.001.fastq" \
#   --output FOLDER=gs://$OUTPUT_B/A549-extract/trimmed_fastq/* \
#   --command 'trim_galore \
#       -a AGATCGGAAGAGCACACGTCTGAAC \
#       -a2 AGATCGGAAGAGCGTCGTGTAGGGA \
#       --quality 30 \
#       --length 40 \
#       --paired \
#       --retain_unpaired \
#       --fastqc \
#       ${R1} \
#       ${R2} \
#       --output_dir $(dirname ${FOLDER})' \
#   --wait



############################## BISMARK ON ONE EXAMPLE

# # Submit job
# dsub \
#   --provider google-v2 \
#   --project $PROJECT_ID \
#   --zones $ZONE_ID \
#   --image $DOCKER_IMAGE \
#   --machine-type n1-highcpu-16 \
#   --disk-size 30 \
#   --logging gs://$OUTPUT_B/logging/ \
#   --input R1="gs://$OUTPUT_B/A549-extract/trimmed_fastq/A549-extract_L01.R1.001_val_1.fq" \
#   --input R2="gs://$OUTPUT_B/A549-extract/trimmed_fastq/A549-extract_L01.R2.001_val_2.fq" \
#   --input-recursive REF_GENOME="gs://$REF_DATA_B/grc37" \
#   --output OUTPUT_DIR="gs://$OUTPUT_B/A549-extract/aligned_per_chard/*" \
#   --command 'bismark_nozip \
#                 -q \
#                 --bowtie2 \
#                 ${REF_GENOME} \
#                 -N 1 \
#                 -1 ${R1} \
#                 -2 ${R2} \
#                 --un \
#                 --score_min L,0,-0.2 \
#                 --bam \
#                 --multicore 3 \
#                 -o $(dirname ${OUTPUT_DIR})' \
#   --wait


####################### Split the chard's BAM per chr

# dsub \
#   --provider google-v2 \
#   --project $PROJECT_ID \
#   --preemptible \
#   --disk-size 100 \
#   --zones $ZONE_ID \
#   --image $DOCKER_IMAGE \
#   --logging gs://$OUTPUT_B/logging/ \
#   --input BAM_FILES="gs://$OUTPUT_B/A549-extract/aligned_per_chard/A549-extract_L01.R1.001_val_1_bismark_bt2_pe.bam" \
#   --output OUTPUT_DIR="gs://$OUTPUT_B/A549-extract/bam_per_chard_and_chr/*" \
#   --script $HOME/GITHUB_REPOS/wgbs-asm/split_bam.sh \
#   --wait



###################### Merge bam per chr

# Submit job
# dsub \
#   --provider google-v2 \
#   --project $PROJECT_ID \
#   --machine-type n1-highmem-8 \
#   --disk-size 30 \
#   --zones $ZONE_ID \
#   --image $DOCKER_IMAGE \
#   --logging gs://$OUTPUT_B/logging/ \
#   --env SAMPLE="A549-extract" \
#   --env CHR="1" \
#   --input BAM_FILES="gs://$OUTPUT_B/A549-extract/bam_per_chard_and_chr/*chr1.bam" \
#   --output OUTPUT_DIR="gs://$OUTPUT_B/A549-extract/bam_per_chr/*" \
#   --script $HOME/GITHUB_REPOS/wgbs-asm/merge_bam.sh \
#   --wait


##################### Net methylation

# dsub \
#   --provider google-v2 \
#   --project $PROJECT_ID \
#   --preemptible \
#   --machine-type n1-highmem-8 \
#   --disk-size 50 \
#   --zones $ZONE_ID \
#   --image $DOCKER_IMAGE \
#   --logging gs://$OUTPUT_B/logging/ \
#   --input BAM="gs://$OUTPUT_B/A549-extract/bam_per_chr/A549-extract_chr1.bam" \
#   --output OUTPUT_DIR="gs://$OUTPUT_B/A549-extract/net_methyl/*" \
#   --command 'bismark_methylation_extractor \
#                   -p \
#                   --no_overlap \
#                   --multicore 3 \
#                   --merge_non_CpG \
#                   --bedGraph \
#                   --counts \
#                   --report \
#                   --buffer_size 48G \
#                   --output $(dirname ${OUTPUT_DIR}) \
#                   ${BAM} \
#                   --ignore 3 \
#                   --ignore_3prime 3 \
#                   --ignore_r2 2 \
#                   --ignore_3prime_r2 2' \
#   --wait



###################### Recal bam

# dsub \
#   --provider google-v2 \
#   --project $PROJECT_ID \
#   --machine-type n1-standard-16 \
#   --disk-size 200 \
#   --zones $ZONE_ID \
#   --image $DOCKER_IMAGE \
#   --logging gs://$OUTPUT_B/logging/ \
#   --env SAMPLE="A549-extract" \
#   --env CHR="1" \
#   --input BAM="gs://$OUTPUT_B/A549-extract/bam_per_chr/A549-extract_chr1.bam" \
#   --input REF_GENOME="gs://$REF_DATA_B/grc37/*" \
#   --input VCF="gs://$REF_DATA_B/dbSNP150_grc37_GATK/no_chr_dbSNP150_GRCh37.vcf" \
#   --output OUTPUT_DIR="gs://$OUTPUT_B/A549-extract/recal_bam_per_chr/*" \
#   --script $HOME/GITHUB_REPOS/wgbs-asm/bam_recalibration.sh \
#   --wait