
#!/bin/bash

# Shortcut to Picardtools
PICARD="/genomics-packages/picard-tools-1.97"

# Create directories
mkdir -p $(dirname "${OUTPUT_DIR}")/$GENOME/variants # For variants
mkdir -p $(dirname "${OUTPUT_DIR}")/$GENOME/ref_genome # For the ref genome

FOLDER=$(dirname "${OUTPUT_DIR}")/$GENOME

# Reference discussion on the ref genome https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use

############################################ Reference genome

if [ $GENOME == "hg19" ] ; then
    GENOME_URL="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz"
else
    # GRCh38
    GENOME_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
fi

# Download zipped ref genome
wget $GENOME_URL

# Unzip file 
gunzip $(basename "$GENOME_URL") 

# Move the file into the output folder
mv $(basename "${GENOME_URL%.gz}") ${FOLDER}/ref_genome/$(basename "${GENOME_URL%.gz}")

# Create the fasta sequence dictionary file (.dict file), required by Bis-SNP
java -jar ${PICARD}/CreateSequenceDictionary.jar \
    R= ${FOLDER}/ref_genome/$(basename "${GENOME_URL%.gz}") \
    O= ${FOLDER}/ref_genome/$(basename "${GENOME_URL%.f*}").dict

# Create a fasta index file, required by bis-snp (*.fai)
samtools faidx ${FOLDER}/ref_genome/$(basename "${GENOME_URL%.gz}")

echo "After dict and fai"
ls -lhR $(dirname "${OUTPUT_DIR}")

#Create bisulfite-converted genome, required by Bismark (Takes ~3 hours)
bismark_genome_preparation --bowtie2 --verbose ${FOLDER}/ref_genome


############################################ Variant database

if [ $GENOME == "hg19" ] ; then
    VARIANTS_URL="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz"
else
    # GRCH38 -- the latest SNP database as of Sept 2019 is db151
    VARIANTS_URL="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz"
fi


# Download zipped variant database
wget $VARIANTS_URL

# Unzip file in the output folder
gunzip $(basename "$VARIANTS_URL")

# Move the file
mv $(basename "${VARIANTS_URL%.gz}") ${FOLDER}/variants/$(basename "${VARIANTS_URL%.gz}")


