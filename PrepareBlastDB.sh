## NOTE: please set your docker disk image size > 120G in Preference>Disk>Disk image size. Or set the blastdb to a attached disk for image.
source /opt/conda/bin/activate jo_ATACseqPipe
cd /blastdb && update_blastdb.pl --decompress nt
