FROM ubuntu:16.04

# To create the image on Cloud Build: gcloud builds submit --timeout "1h" --tag gcr.io/hackensack-tyco/onuchic-perl

# Update package list
RUN apt-get update -y 

# Bundle of useful packages (Make is necessary to install cpanminus, etc.)
RUN apt-get install build-essential -y

# wget is not included in the essentials
RUN apt-get install wget -y

# Install Perl modules (e.g cpan)
RUN apt-get install perl-modules -y

# Install cpanm
RUN cpan App::cpanminus

# Install samtools 0.1.19 (newer versions are NOT compatible)
RUN wget 'https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2'
RUN tar -xvf samtools-0.1.19.tar.bz2
RUN apt-get --yes install zlib1g-dev libncurses5-dev libncursesw5-dev 
RUN cd samtools-0.1.19 && make
ENV PATH="samtools-0.1.19:${PATH}"

# Install samtools for Perl
RUN cpanm Bio::DB::Sam

##################### INSTALLATION END #####################

# Default job to execute, here: display the date.
CMD ["date"]
