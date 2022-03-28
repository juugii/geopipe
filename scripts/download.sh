#! /usr/bin/env sh

# Create a list of SRR runs from a GDS id
# Jules GILET <jules.gilet@inserm.fr>

# GDS id as arg1
# output file as arg2

echo "Fetching SRRs..."

for I in $(esearch -db gds -query ${1} | efetch -format native | grep -e "SRX"| cut -d "=" -f 2); do
	esearch -db sra -query ${I} | efetch -format runinfo | grep -e "SRR" | cut -d "," -f 1
done > ${2}

echo "Created SRR list file."

exit 0
