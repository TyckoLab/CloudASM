FROM bioconductor/release_base2:R3.6.1_Bioc3.9
# It is a Debian machine

# To create the image on Cloud Build: gcloud builds submit --timeout "1h" --tag gcr.io/hackensack-tyco/r-genomics

# Update package list
RUN apt-get update -y 

# R Package within Bioconductor
RUN R -e "BiocManager::install('GenomicRanges')"

##################### INSTALLATION END #####################

# Default job to execute, here: display the date.
CMD ["date"]
