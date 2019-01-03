From:r-base:latest
Bootstrap:docker

%labels
    MAINTAINER Jianhong Ou <jianhong.ou@duke.edu>
    DESCRIPTION Singularity image containing all requirements for the ATACseqQCnextflow
    VERSION 0.0.3

%files
    condaEnv.yml /
    src /tmp/
    blastDB/ncbirc /etc/.ncbirc

%post
  apt-get update --fix-missing
  apt-get install --yes wget git bzip2 ca-certificates curl unzip gdebi-core rsync libssl-dev libcurl4-openssl-dev libgsl-dev zlib1g-dev ttf-dejavu g++ libxml2 libxml2-dev ghostscript

  cd /tmp/src && g++ -o /usr/local/bin/fq2sc fq2sc.cpp -lz
  cp s2c_2_fasta.pl /usr/local/bin/ && chmod +x /usr/local/bin/s2c_2_fasta.pl

  mkdir -p /blastdb && mkdir -p /igenome
  echo 'install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)' | R --vanilla
  echo 'BiocManager::install("ChIPpeakAnno")' | R --vanilla
  echo 'BiocManager::install(c("GenomicScores", "randomForest", "motifStack", "preseqR"))' | R --vanilla
  cd /tmp/ && git clone https://github.com/jianhong/ATACseqQC.git && R CMD INSTALL ATACseqQC

  /opt/conda/bin/conda env create -f /condaEnv.yml
  /opt/conda/bin/conda clean -a
  /opt/conda/bin/activate jo_ATACseqPipe
