#!/bin/bash
usage="$(basename "$0") [-h] [-spfnetcugd] -- RNAseq_pipeline

where:
    -h  show this help text
    -s  species; default danRer10
    -p  postfix of fastq files; default fastq.gz
    -f  folder where save fastq files; default fastq
    -n  surname for the bash files; default atacseq
    -e  single end (se) or paired end (pe); default se
    -c  number of threads; default 1
    -t  temp folder; default /tmp
    -u  root of reference files; default iGenomes. The subfolder should be species/Sequence/...
    -g  genomeSize for MACS2 (--gsize)
    -d  dry-run
	
example:
	$(basename "$0") -s danrer10 -p fastq.gz
	
Note: 
	fastq.gz file must be in a subfolder for example: fastq/inj or fastq/uninj.
	Each fastq file should be one biological replicate.
	If you have multiple fastq files for one replicate, please cat then together.
"


prefix=bowtie2
postfix=fastq.gz
folder=fastq
species=danRer10
surname=atacseq
singleEnd=se
cups=1
tmpfolder="/tmp"
UCSC="iGenomes"
run="yes"
genome=0

decleare -A genomeSize=(["hg38"]="2.7e9", ["mm10"]="1.87e9", ["ce6"]="9e7", ["dm3"]="1.2e8", ["danRer10"]="1.5e9")

while getopts ':hs:p:f:n:t:d' option; do
	case "$option" in
		h) echo "$usage"
		   exit
		   ;;
		s) species=$OPTARG
		   ;;
		p) postfix=$OPTARG
		   ;;
		f) folder=$OPTARG
		   ;;
		n) surname=$OPTARG
		   ;;
		e) singleEnd=$OPTARG
		   ;;
		c) cpus=$OPTARG
		   ;;
		u) UCSC=$OPTARG
		   ;;
		g) genome=$OPTARG
		   ;;
		d) run="no"
		   ;;
		:) printf "missing argument for -%s\n" "$OPTARG" >&2
       	   echo "$usage" >&2
           exit 1
           ;;
		\?) printf "illegal option: -%s\n" "$OPTARG" >&2
           echo "$usage" >&2
           exit 1
           ;;
	esac
done
shift $((OPTIND - 1))

if [ "$singleEnd" != "se" ] && [ "$singleEnd" != "pe" ]; then
	printf "-e argument must be se or pe" >&2
	echo "$usage" >&2
	exit
fi

if [ "$singleEnd" == "pe" ]; then
	printf "only support se now" >&2
	echo "$usage" >&2
	exit
fi

if [ "$genome" -eq 0 ]; then
  genome=${genomeSize[$species]}
fi

files=(`find $folder/ -iname "*.$postfix"`)
group=()
tag=()
filenames=()
for i in ${files[@]}
do
	this=`basename $i`
	this=`echo $this | sed -e "s/.$postfix//"`
	filenames+=($this)
	this=`dirname $i`
	this=`basename $this`
	group+=($this)
done

labels=($(echo "${group[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

declare -A lab
for i in ${labels[@]}
do
	lab+=([$i]=0)
done

for i in ${group[@]}
do
	lab[$i]=$((${lab[$i]} + 1))
	tag+=($i.rep${lab[$i]})
done

mkdir -p log
mkdir -p fastqc
mkdir -p fastq.trimmed
mkdir -p sam
mkdir -p bam
mkdir -p macs2
mkdir -p bigwig
mkdir -p taxReport

numOfLine=100000
mismatch=2
numOfL=$(($numOfLine * 4))

for i in ${filenames[@]}
do
  ## fastqc
  mkdir -p fastqc/${group[$i]}
  fastqc -o fastqc/${group[$i]} -t $cpus \
         $folder/${group[$i]}/${filenames[$i]}.$postfix
  
  ## trim adapter
  mkdir -p fastq.trimmed/${group[$i]}
  trim_galore -q 15 --fastqc -o fastq.trimmed/${group[$i]} $folder/${group[$i]}/${filenames[$i]}.$postfix
  
  ## taxReport
  zcat fastq.trimmed/${group[$i]}/${filenames[$i]}_trimmed.fq.gz | head -n $numOfL > taxReport/${filenames[$i]}.fq
  fq2sc -i taxReport/${filenames[$i]}.fq -o taxReport/${filenames[$i]}.txt
  s2c_2_fasta.pl --in taxReport/${filenames[$i]}.txt --out taxReport/${filenames[$i]}.fa --notfilterN
  blastn -num_threads $cpus -db nt -max_target_seqs 1 \
  -outfmt '6 qaccver saccver staxid sblastname ssciname scomname pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
  -query taxReport/${filenames[$i]}.fa -out taxReport/${filenames[$i]}.blastn.txt
  cat <<ET | R --vanilla
x <- read.delim("taxReport/${filenames[$i]}.blastn.txt", header=FALSE, stringsAsFactor=FALSE)
colnames(x) <- c("qaccver", "saccver", "staxid", "sblastname", "ssciname", "scomname", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
x <- x[x\\\$mismatch<=$mismatch, ]
x <- x[order(x\\\$qaccver, -x\\\$pident, x\\\$evalue), ]
x <- x[!duplicated(x\\\$qaccver), ]
sc <- as.numeric(sub("SEQ_.*?_x", "", x\\\$qaccver))
species <- paste0(x\\\$ssciname, " (", x\\\$scomname, ")", " [", x\\\$sblastname, "]")
w <- rowsum(sc, species)
w <- w[order(-w[, 1]), ]
p <- w/$numOfLine
write.table(p, 'taxReport/${filenames[$i]}.percentage.txt', quote=FALSE, col.names=FALSE, sep='\t')
f <- dir("taxReport", "percentage.txt", full.names = TRUE)
d <- lapply(f, read.delim, nrows=3, header=FALSE)
e <- sapply(d, function(.ele){
  if(nrow(.ele)>1){
    .ele[1, 2]/sum(.ele[-1, 2])
  }else{
    100
  }
})
e <- cut(e, breaks=c(0, 3, 5, 10, Inf), labels=c("contaminated", "warning", "minor", "clean"))
e <- rep(e, each=3)
fq <- sub(".$postfix.percentage.txt", "", basename(f))
fq <- rep(fq, each=3)
d <- do.call(rbind, d)
d <- cbind(d, fastq=fq, cross.species.contamination=e)
colnames(d) <- c("species", "percentage", "fastq", "cross.species.contamination")
write.csv(d, "taxReport.csv", row.names=FALSE)
ET
  rm taxReport/${filenames[$i]}.fa
  rm taxReport/${filenames[$i]}.txt
  rm taxReport/${filenames[$i]}.fastq
  
  ## mapping by bowtie2
  bowtie2 -x $UCSC/$species/Sequence/Bowtie2Index/genome \
        -U fastq.trimmed/${group[\$i]}/${filenames[$i]}_trimmed.fq.gz \
        -S sam/$prefix.$species.${tag[$i]}.sam \
        --very-sensitive \
        -p $cpus \
        > sam/$prefix.$species.${tag[$i]}.log
  
  ## sam to bam
  samtools view -bhS -q 30 sam/$prefix.$species.${tag[$i]}.sam > bam/$prefix.$species.${tag[$i]}.bam
  rm sam/$prefix.$species.${tag[$i]}.sam
  samtools sort bam/$prefix.$species.${tag[$i]}.bam -o bam/$prefix.$species.${tag[$i]}.srt.bam
  samtools index bam/$prefix.$species.${tag[\$i]}.srt.bam
  rm bam/$prefix.$species.${tag[$i]}.bam
  
  ## remove PCR duplicates
  picard MarkDuplicates INPUT=bam/$prefix.$species.${tag[$i]}.srt.bam \
	OUTPUT=bam/markDup/$prefix.$species.${tag[$i]}.markDup.bam \
	METRICS_FILE=bam/$prefix.$species.${tag[$i]}.srt.fil.picard_info.txt \
	REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
  
  ## call peaks
  mkdir -p macs2/${tag[$i]}
  macs2 callpeak -t bam/$prefix.$species.${tag[$i]}.markDup.bam \
                 -f BAM -g $genome  -n $prefix.$species.${tag[$i]} \
                 --outdir macs2/${tag[$i]} -q 0.05 \
                 --nomodel --shift 37 --extsize 73 -B


  # generate bigwig files
  curr=$prefix.$species.${tag[$i]}.markDup
  samtools view -H bam/$curr.bam | grep @SQ | awk -F $'\t' 'BEGIN {OFS=FS} {gsub("[SL]N:", ""); print $2,$3}' > $curr.$i.ChromInfo.txt
  sort -k1,1 -k2,2n macs2/${tag[$i]}/$prefix.$species.${tag[$i]}_treat_pileup.bdg macs2/${tag[$i]}/$prefix.$species.${tag[$i]}_treat_pileup.srt.bdg
  bedGraphToBigWig macs2/${tag[$i]}/$prefix.$species.${tag[$i]}_treat_pileup.srt.bdg $curr.$i.ChromInfo.txt bigwig/$prefix.$species.${tag[$i]}.bw
  rm $curr.$i.ChromInfo.txt
  rm macs2/${tag[$i]}/$prefix.$species.${tag[$i]}_treat_pileup.srt.bdg
done

