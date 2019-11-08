FROM bioperl/bioperl

# To create the image on Cloud Build: gcloud builds submit --timeout "1h" --tag gcr.io/hackensack-tyco/perl-samtools

# Update package list
RUN apt-get update -y 

# wget is not included in the essentials
RUN apt-get install wget -y

# Install samtools 0.1.19 (newer versions are NOT compatible)
RUN apt-get --yes install gcc libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev
RUN wget "https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2"
RUN tar xjf samtools-0.1.19.tar.bz2
WORKDIR /samtools-0.1.19
RUN make CFLAGS=-fPIC
ENV SAMTOOLS="/samtools-0.1.19/"
ENV PATH="/samtools-0.1.19/:${PATH}"
WORKDIR /

# Install samtools for Perl
RUN cpanm Bio::DB::Sam

##################### INSTALLATION END #####################

# Default job to execute, here: display the date.
CMD ["date"]



