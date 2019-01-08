From:r-base:latest
Bootstrap:docker

%labels
    MAINTAINER Jianhong Ou <jianhong.ou@duke.edu>
    DESCRIPTION Singularity image containing all requirements for the ATACseqQCnextflow
    VERSION 0.0.3

%post
  apt-get update --fix-missing
  apt-get install --yes wget git bzip2 ca-certificates curl unzip gdebi-core rsync libssl-dev libcurl4-openssl-dev libgsl-dev zlib1g-dev ttf-dejavu g++ libxml2 libxml2-dev ghostscript
  
  tmpfold=`pwd`
  git clone https://github.com/jianhong/ATACseqPipe.git
  cd $tmpfold/ATACseqPipe/src && g++ -o /usr/local/bin/fq2sc fq2sc.cpp -lz
  cp s2c_2_fasta.pl /usr/local/bin/ && chmod +x /usr/local/bin/s2c_2_fasta.pl
  cp ../blastDB/ncbirc /etc/.ncbirc

  mkdir -p /blastdb && mkdir -p /igenome
  echo 'install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)' | R --vanilla
  echo 'BiocManager::install("ChIPpeakAnno")' | R --vanilla
  echo 'BiocManager::install(c("GenomicScores", "randomForest", "motifStack", "preseqR"))' | R --vanilla
  cd $tmpfold && git clone https://github.com/jianhong/ATACseqQC.git && R CMD INSTALL ATACseqQC


## add conda
  wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
    

## test conda
  /opt/conda/bin/conda update -y conda

#  /opt/conda/bin/conda env create -f $tmpfold/ATACseqPipe/condaEnv.yml
#  /opt/conda/bin/conda clean -a
#  /opt/conda/bin/activate jo_ATACseqPipe

  /opt/conda/bin/conda install -y -c bioconda fastqc python trim-galore blast bowtie2 \
  bwa samtools picard macs2 bedtools deeptools ucsc-bedgraphtobigwig deeptools
## fix the issue of samtools cannot find libcrypto.so.1.0.0
  ln -s /opt/conda/lib/libcrypto.so.1.1 /opt/conda/lib/libcrypto.so.1.0.0
  ln -s /lib/x86_64-linux-gnu/libssl.so.1.1 /lib/x86_64-linux-gnu/libssl.so.1.0.0

  /bin/bash -c "echo 'BiocManager::install(\"MotifDb\")' | R --vanilla"
