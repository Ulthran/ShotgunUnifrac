#!/bin/bash
# Downloads genome from NCBI for a given genome id
# @param OUT is the output path for the fasta file in the form
#        dir/ID_restOfPath.fasta

OUT=$1
ID=${OUT:5:-23} # Relies on there being 5 leading characters and 23 following around the id
echo $ID
url_prefix=https://ftp.ncbi.nlm.nih.gov/genomes/all
url=$url_prefix/${ID:0:3}/${ID:4:3}/${ID:7:3}/${ID:10:3}/$ID
file_name=${ID}_cds_from_genomic.fna.gz
cds_url=$url/$file_name
echo $cds_url

cd ncbi/
[ -e ${file_name::-7}.fasta ] || wget $cds_url
gzip -d $file_name
mv ${file_name::-3} ${file_name::-7}.fasta
cd ..