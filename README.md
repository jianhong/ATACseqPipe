# ATACseqPipe
ATAC-seq analysis Pipeline


## Install

### Docker

Please set your docker disk image size > 120G in Preference>Disk>Disk image size. Or set the blastdb to a attached disk for image.

```{bash}
cd ~
docker pull jianhong/atacseqpipe:latest
mkdir tmp4atacseqpipe
docker run -it --rm -v ${PWD}/tmp4atacseqpipe:/volume/data \
       jianhong/atacseqpipe:latest bash

## In the docker
cd /volume/data
git clone https://github.com/jianhong/ATACseqPipe.git
cd ATACseqPipe
bash PrepareBlastDB.sh
bash PrepareGenome.sh -s danRer10
```

## Run

### Docker

Please copy your data into tmp4atacseqpipe folder in following format:

```
fastq-|
      |--condition1-|
      |             |--rep1.fastq.gz
      |             |--rep2.fastq.gz
      |             |-- ...
      |             |--repN.fastq.gz
      |--condition2-|
      |             |--rep1.fastq.gz
      |             |--rep2.fastq.gz
      |             |-- ...
      |             |--repN.fastq.gz
      |...
```

And in the docker:

```{bash}
source /opt/conda/bin/activate jo_ATACseqPipe
bash ATACseq.sh
```

