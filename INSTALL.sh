#!/bin/bash
usage="$(basename "$0") [-h] 

where:
  -h  show this help text
  -p  prefix for conda. Default is /opt/conda

example:
  $(basename "$0") -p "$HOME"/conda 
"
prefix="/opt/conda"
while getopts ':hp:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    p) prefix=$OPTARG
       ;;
  esac
done

## Insatll conda
wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    /bin/bash miniconda.sh -b -p $prefix && \
    rm miniconda.sh && \
    $prefix/bin/conda clean -tipsy && \
    echo ". $prefix/etc/profile.d/conda.sh" >> ~/.bashrc

## Install the dependencies:
$prefix/bin/conda env create -f condaEnv.yml

## activate the environment:
source $prefix/bin/activate jo_ATACseqPipe

## install Bioconductor
echo "install.packages(\"BiocManager\", repos='https://cloud.r-project.org', quiet=TRUE)" | R --vanilla
#echo "BiocManager::install(c(\"TxDb.Hsapiens.UCSC.hg38.knownGene\", \"org.Hs.eg.db\", \"TxDb.Drerio.UCSC.danRer10.refGene\", \"org.Dr.eg.db\", \"WriteXLS\", \"ggrepel\"), suppressUpdates=TRUE, ask=FALSE)" | R --vanilla
#echo "BiocManager::install(c(\"ChIPpeakAnno\", \"trackViewer\", \"motifStack\", \"ATACseqQC\", \"GeneNetworkBuilder\", \"DESeq2\", \"csaw\", \"DiffBind\"), suppressUpdates=TRUE, ask=FALSE)" | R --vanilla

## download the blastdb
mkdir -p blastdb
cd blastdb && update_blastdb.pl --decompress nt

## prepare for taxReport
wget -O $prefix/bin/s2c_2_fasta.pl https://raw.githubusercontent.com/jianhong/ATACseqPipe/master/src/s2c_2_fasta.pl
wget https://raw.githubusercontent.com/jianhong/ATACseqPipe/master/src/fq2sc.cpp
wget https://raw.githubusercontent.com/jianhong/ATACseqPipe/master/src/fq2sc.h
wget https://raw.githubusercontent.com/jianhong/ATACseqPipe/master/src/kseq.h
g++ -o $prefix/bin/fq2sc fq2sc.cpp -lz

## download the reference files
mkdir -p igenome
echo "Please download reference files from https://support.illumina.com/sequencing/sequencing_software/igenome.html into igenome folder"
