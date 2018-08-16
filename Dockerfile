#################################################################
# Dockerfile to build bowtie2, MACS2, samtools, 
# picard-tools, fastQC, bedtools, cutadapt, R, blast
# images
# Based on Ubuntu
#  $ cd ATACseqPipe.docker
#  $ VERSION=0.0.1
#  $ docker build -t jianhong/ATACseqPipe:$VERSION .
#  $ docker images jianhong/ATACseqPipe:$VERSION
#  $ docker push jianhong/ATACseqPipe:$VERSION
#  $ docker tag jianhong/ATACseqPipe:$VERSION jianhong/ATACseqPipe:latest
#  $ docker push jianhong/ATACseqPipe:latest
#  $ cd ~
#  $ docker pull jianhong/ATACseqPipe:latest
#  $ mkdir tmp4ATACseqPipe
#  $ docker run -it --rm -v ${PWD}/tmp4ATACseqPipe:/volume/data \
#  $       jianhong/ATACseqPipe:latest bash
##################################################################
# Set the base image to Ubuntu
FROM ubuntu:latest

# File/Author / Maintainer
MAINTAINER Jianhong Ou <jianhong.ou@duke.edu>

# envirenment
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH $PATH:/opt/conda/bin

# Update the repository sources list, install wget, unzip, curl, git
RUN \
  apt-get update --fix-missing && \
  apt-get install --yes wget git bzip2 ca-certificates curl unzip gdebi-core git rsync libssl-dev libcurl4-openssl-dev libgsl-dev zlib1g-dev g++ && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*
  
## add conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
    

## test conda
RUN /opt/conda/bin/conda update -y conda

## git jianhong/ATACseqPipe and install the dependencies:
RUN git clone https://github.com/jianhong/ATACseqPipe.git && \
    cd ATACseqPipe && /opt/conda/bin/conda env create -f condaEnv.yml && \
    source /opt/conda/bin/activate jo_ATACseqPipe && \
    cd src && g++ -o /usr/local/bin/fq2sc fq2sc.cpp -lz && \
    cp s2c_2_fasta.pl /usr/local/bin/ && chmod +x /usr/local/bin/s2c_2_fasta.pl && \
    cp ../blastDB/ncbirc /etc/.ncbirc && \
    cd / && echo "install.packages(\"BiocManager\", repos='https://cloud.r-project.org', quiet=TRUE)" | R --vanilla

## make directory
RUN mkdir -p /blastdb && mkdir -p /igenome


