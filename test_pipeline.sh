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




######## TEST GENOMIC OVERLAP WITH BIG QUERY



# Import the VCF file in Big Query without the VCF header
bq --location=US load \
               --replace=false \
               --source_format=CSV \
               --field_delimiter "\t" \
               ${DATASET_ID}.hg19_cpg_pos \
               gs://$REF_DATA_B/hg19_CpG_pos.bed \
               chr:STRING,inf:INT64,sup:INT64

# Creates duplicate rows with snp_id

WITH
  variants AS (
    SELECT
      chr AS chr_variants,
      inf AS inf_variants,
      sup AS sup_variants
    FROM
      `hackensack-tyco.wgbs_asm.hg19_cpg_pos`
  )
  SELECT DISTINCT
    *
  FROM
     `hackensack-tyco.wgbs_asm.gm12878_21_500bp` AS intervals
  INNER JOIN
    variants ON
      variants.chr_variants = chr
      AND variants.inf_variants >= inf_500
      AND variants.sup_variants <= sup_500
  GROUP BY 
      snp_id


# Working

WITH
  variants AS (
    SELECT
      chr AS chr_variants,
      inf AS inf_variants,
      sup AS sup_variants
    FROM
      `hackensack-tyco.wgbs_asm.hg19_cpg_chr21_extract`
  )
  SELECT DISTINCT
    snp_id,
    chr,
    pos,
    coord_2000,
    ref,
    alt,
    cov
  FROM
     `hackensack-tyco.wgbs_asm.gm12878_21_500bp`
  INNER JOIN
    variants ON
      variants.chr_variants = chr
      AND variants.inf_variants >= inf_500
      AND variants.sup_variants <= sup_500



#### WORKING

WITH
  variants AS (
    SELECT
      chr AS chr_variants,
      inf AS inf_variants,
      sup AS sup_variants
    FROM
      `hackensack-tyco.wgbs_asm.hg19_cpg_chr21_extract`
  )
  SELECT DISTINCT
    snp_id,
    chr,
    pos,
    coord_2000,
    ref,
    alt,
    cov
  FROM
     `hackensack-tyco.wgbs_asm.gm12878_21_500bp`
  INNER JOIN
    variants ON
      variants.chr_variants = chr
      AND variants.inf_variants >= inf_500
      AND variants.sup_variants <= sup_500



######################### Big query working but not reliable.

bq query \
    --use_legacy_sql=false \
    --destination_table ${PROJECT_ID}:${DATASET_ID}.${SAMPLE}_chr${CHR}_500bp_wCpG \
    --replace=true \
    "WITH
      win_2000 AS (
        SELECT 
          chr,
          snp_id,
          SAFE_CAST(pos AS INT64) - 250 AS inf_500, 
          SAFE_CAST(pos AS INT64) + 250 AS sup_500,
          SAFE_CAST(pos AS INT64) - 1000 AS inf_2000,
          SAFE_CAST(pos AS INT64) + 1000 AS sup_2000,
          ref,
          alt,
          SAFE_CAST(REGEXP_EXTRACT(info,'DP=(.+);M') as INT64) as cov,
          SAFE_CAST(pos as INT64) as pos
        FROM 
          ${DATASET_ID}.${SAMPLE}_chr${CHR}_vcf    
        ),
      -- we convert the 2000bp window boundaries into string
      win_500 AS (
        SELECT 
          chr, 
          inf_500,
          sup_500,
          snp_id,
          ref,
          alt,
          cov,
          pos,
          SAFE_CAST(inf_2000 AS STRING) as inf_2000,
          SAFE_CAST(sup_2000 AS STRING) as sup_2000
        FROM win_2000),
      -- create a coordinate of a 2000 bp window (to be used later)
      variants AS (
        SELECT 
          chr,
          inf_500,
          sup_500,
          snp_id,
          ref,
          alt,
          cov,
          pos,
          CONCAT(chr,':', inf_2000, '-', sup_2000) as coord_2000
        FROM win_500
        -- we ask that snp ID are rs* and that ref and alt nucleotides have only
        WHERE 
          snp_id LIKE '%rs%'
          AND cov >= 10
          AND BYTE_LENGTH(ref) = 1
          AND BYTE_LENGTH(alt) = 1),
      cpg_pos AS (
        SELECT
          chr AS chr_cpg_pos,
          inf AS inf_cpg_pos,
          sup AS sup_cpg_pos
        FROM
          ${DATASET_ID}.hg19_cpg_pos)
  -- We make sure that there is a CpG in the 500bp window of the SNP
  SELECT DISTINCT
    snp_id,
    chr,
    pos,
    coord_2000,
    ref,
    alt,
    cov
  FROM
     variants
  INNER JOIN
    cpg_pos ON
      cpg_pos.chr_cpg_pos = chr
      AND cpg_pos.inf_cpg_pos >= inf_500
      AND cpg_pos.sup_cpg_pos <= sup_500
  "


############## Load CpG OB and OT of chr 22

bq --location=US load \
               --replace \
               --source_format=CSV \
               --field_delimiter "\t" \
               --skip_leading_rows 117 \
               ${DATASET_ID}.${SAMPLE}_chr${CHR}_vcf \
               gs://$BUCKET/$SAMPLE/variants_per_chr/${SAMPLE}_chr${CHR}.vcf \
               chr:STRING,pos:STRING,snp_id:STRING,ref:STRING,alt:STRING,qual:FLOAT,filter:STRING,info:STRING,format:STRING,data:STRING
