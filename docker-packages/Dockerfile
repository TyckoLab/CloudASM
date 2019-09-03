FROM google/cloud-sdk:255.0.0

## To create the image on Cloud Build: gcloud builds submit --timeout "1h" --tag gcr.io/hackensack-tyco/wgbs-asm

RUN apt-get update -y 
# We do NOT upgrade: apt-get upgrade -y
RUN apt-get install build-essential -y

# Convenient text editor
RUN apt-get install nano

################## BEGIN INSTALLATION OF GENOMIC PACKAGES ######################

WORKDIR /root

ARG PACKAGES="/genomics-packages"
RUN mkdir -p $PACKAGES

# Bowtie (used by Samtools)
RUN gsutil -m cp -r gs://genomic-packages/bowtie2-2.1.0 $PACKAGES
RUN cd $PACKAGES/bowtie2-2.1.0 && make
ENV PATH="${PACKAGES}/bowtie2-2.1.0:${PATH}"

# Samtools (manipulation of BAM and SAM files)
RUN gsutil -m cp -r gs://genomic-packages/samtools-0.1.19 $PACKAGES
RUN apt-get --yes install zlib1g-dev libncurses5-dev libncursesw5-dev 
RUN cd $PACKAGES/samtools-0.1.19 && make
ENV PATH="${PACKAGES}/samtools-0.1.19:${PATH}"

# Bedtools (to manipulate bedgraph files)
RUN gsutil -m cp -r gs://genomic-packages/bedtools2 $PACKAGES
ENV PATH="${PACKAGES}/bedtools2/bin:${PATH}"

# Bismark (Aligner)
RUN gsutil -m cp -r gs://genomic-packages/Bismark_v0.19.0 $PACKAGES
ENV PATH="${PACKAGES}/Bismark_v0.19.0:${PATH}"

# Cutadapt
RUN gsutil -m cp -r gs://genomic-packages/cutadapt-1.4.2 $PACKAGES
ENV PATH="${PACKAGES}/cutadapt-1.4.2/bin:${PATH}"

#TrimGalore
RUN gsutil -m cp -r gs://genomic-packages/TrimGalore-0.4.5 $PACKAGES
ENV PATH="${PACKAGES}/TrimGalore-0.4.5:${PATH}"

# PicardTools
RUN gsutil -m cp -r gs://genomic-packages/picard-tools-1.97/picard-tools-1.97 $PACKAGES
ENV PATH="${PACKAGES}/picard-tools-1.97:${PATH}"

# Java (Bis-SNP requires a specific version to run properly)
RUN gsutil -m cp -r gs://genomic-packages/jre1.7.0_25 $PACKAGES
ENV PATH="${PACKAGES}/jre1.7.0_25:${PATH}"
ENV PATH="${PACKAGES}/jre1.7.0_25/bin:${PATH}"

# Bis-SNP (variant call in WGBS data)
RUN gsutil -m cp -r gs://genomic-packages/BisSNP-0.82.2 $PACKAGES
ENV PATH="${PACKAGES}/BisSNP-0.82.2:${PATH}"

# FastQC
RUN gsutil -m cp -r gs://genomic-packages/FastQC $PACKAGES
ENV PATH="${PACKAGES}/FastQC:${PATH}"

# Gives rights to call the packages
RUN chmod -R 777 $PACKAGES

####################  INSTALL R  ######################

# Install R 
RUN apt-get install r-base r-base-dev -y

# Multi-thread for linear operations in R
RUN apt-get install libatlas3-base -y

# R packages used by scripts
RUN su - -c "R -e \"install.packages('data.table', repos='http://cran.rstudio.com/')\""

##################### INSTALLATION END #####################

# Default job to execute, here: display the date.
CMD ["date"]

# Set default container command (this will avoid Docker exit code 2 when run locally.)
#ENTRYPOINT ["/bin/sh -c"]