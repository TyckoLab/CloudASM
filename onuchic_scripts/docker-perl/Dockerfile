FROM ubuntu:latest

## To create the image on Cloud Build: gcloud builds submit --timeout "1h" --tag gcr.io/hackensack-tyco/wgbs-asm

RUN apt-get update -y 
RUN apt-get install build-essential -y

################## BEGIN INSTALLATION OF GENOMIC PACKAGES ######################

WORKDIR /root


##################### INSTALLATION END #####################

# Default job to execute, here: display the date.
CMD ["date"]

# Set default container command (this will avoid Docker exit code 2 when run locally.)
#ENTRYPOINT ["/bin/sh -c"]