#################################################################
# Dockerfile to build bowtie2, MACS2, samtools, 
# picard-tools, fastQC, bedtools, cutadapt, R, blast
# images
# Based on r-base:latest
#  $ cd ATACseqPipe.docker
#  $ VERSION=0.0.3
#  $ docker build --no-cache -t jianhong/atacseqpipe:$VERSION .
#  $ docker images jianhong/atacseqpipe:$VERSION
#  $ docker push jianhong/atacseqpipe:$VERSION
#  $ docker tag jianhong/atacseqpipe:$VERSION jianhong/atacseqpipe:latest
#  $ docker push jianhong/atacseqpipe:latest
#  $ cd ~
#  $ docker pull jianhong/atacseqpipe:latest
#  $ mkdir tmp4atacseqpipe
#  $ docker run -it --rm -v ${PWD}/tmp4atacseqpipe:/volume/data \
#  $       jianhong/atacseqpipe:latest bash
##################################################################
# Set the base image to Ubuntu
FROM r-base:latest

# File/Author / Maintainer
MAINTAINER Jianhong Ou <jianhong.ou@duke.edu>

# envirenment
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH $PATH:/opt/conda/bin

# Update the repository sources list, install wget, unzip, curl, git
RUN \
  apt-get update --fix-missing && \
  apt-get install --yes wget git bzip2 ca-certificates curl unzip gdebi-core rsync libssl-dev libcurl4-openssl-dev libgsl-dev zlib1g-dev ttf-dejavu g++ libxml2 libxml2-dev && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

# install fq2sc and s2c_2_fasta.pl
RUN cd /tmp/ && git clone https://github.com/jianhong/ATACseqPipe.git && \
    cd /tmp/ATACseqPipe/src && g++ -o /usr/local/bin/fq2sc fq2sc.cpp -lz && \
    cp s2c_2_fasta.pl /usr/local/bin/ && chmod +x /usr/local/bin/s2c_2_fasta.pl && \
    cp ../blastDB/ncbirc /etc/.ncbirc

## make directory
RUN mkdir -p /blastdb && mkdir -p /igenome

## install ATACseqQC
##RUN /bin/bash -c "source /opt/conda/bin/activate jo_ATACseqPipe && echo 'BiocManager::install(\"ATACseqQC\", version=\"3.8\")' | R --vanilla"

RUN /bin/bash -c "echo 'install.packages(\"BiocManager\", repos=\"https://cloud.r-project.org\", quiet=TRUE)' | R --vanilla"
RUN /bin/bash -c "echo 'BiocManager::install(\"ChIPpeakAnno\")' | R --vanilla"
RUN /bin/bash -c "echo 'BiocManager::install(c(\"GenomicScores\", \"randomForest\", \"motifStack\", \"preseqR\"))' | R --vanilla"
RUN cd /tmp/ && git clone https://github.com/jianhong/ATACseqQC.git && \
  R CMD INSTALL ATACseqQC
#RUN /bin/bash -c "echo 'BiocManager::install(\"jianhong/ATACseqQC\")' | R --vanilla"

## add conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
    

## test conda
RUN /opt/conda/bin/conda update -y conda && \
  /opt/conda/bin/conda install -y -c bioconda fastqc python trim-galore blast bowtie2 \
  bwa samtools picard macs2 bedtools deeptools ucsc-bedgraphtobigwig deeptools
## fix the issue of samtools cannot find libcrypto.so.1.0.0
RUN ln -s /opt/conda/lib/libcrypto.so.1.1 /opt/conda/lib/libcrypto.so.1.0.0
RUN ln -s /lib/x86_64-linux-gnu/libssl.so.1.1 /lib/x86_64-linux-gnu/libssl.so.1.0.0

RUN /bin/bash -c "echo 'BiocManager::install(\"MotifDb\")' | R --vanilla"