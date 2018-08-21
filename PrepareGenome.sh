#!/bin/bash
usage="$(basename "$0") [-h] [-s]

where:
    -h  show this help text
    -s  species; default danRer10
    
example:
	$(basename "$0") -s hg38
"

species=danRer10

while getopts ':hs:' option; do
	case "$option" in
		h) echo "$usage"
		   exit
		   ;;
		s) species=$OPTARG
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

declare -A sciencename=(["hg38"]="Homo_sapiens", ["hg19"]="Homo_sapiens", \
                         ["mm10"]="Mus_musculus", ["ce6"]="Caenorhabditis_elegans", \
                         ["dm3"]="Drosophila_melanogaster", ["danRer10"]="Danio_rerio")

cd /igenome && \
 wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/${sciencename[$species]}/UCSC/${species}/${sciencename[$species]}_UCSC_${species}.tar.gz && \
 tar -xzf ${sciencename[$species]}_UCSC_${species}.tar.gz && rm ${sciencename[$species]}_UCSC_${species}.tar.gz
 
 