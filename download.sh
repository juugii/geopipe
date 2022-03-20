#! /bin/env sh

Jules GILET <jules.gilet@inserm.fr>


GDS="GSE135194"

# query a GEO ID and retrieve SRX
SRX=$(esearch -db gds -query ${GDS} | efetch -format native | grep SRX | cut -d "=" -f 2)

# translate SRX to SRR
SRR=$(for A in ${SRX[@]}; do esearch -db sra -query ${A} | efetch -format runinfo | grep -e "SRR" | cut -d "," -f 1; done)

# sra download and seqfile dump
for SRA in ${SRR[@]}; do
	prefetch --max-size 100GB ${SRA} && vdb-validate ${SRA}
done && fastq-dump --split-files */*.sra

# --gzip is deprecated
for FILE in *.fastq; do
	# file 1 is index
	# read 1 is file 2
	mv ${FILE} ${FILE//_2/_S1_L001_R1_001}
	# read 2 is file 3
	mv ${FILE} ${FILE//_3/_S1_L001_R2_001}
done

# compression (dev dataset is 330 Gio..)
for FILE in *001.fastq; do 
	gzip ${FILE}
done

exit 0
