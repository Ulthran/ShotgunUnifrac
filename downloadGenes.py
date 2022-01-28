#!~/mambaforge/envs/shotguntree_env/bin/ python3
# Downloads genome from NCBI for a given genome id
# @param out is the output file continaing the accession of the desired genome for download
#        taking the form ncbi/{ID}_cds_from_genomic.fasta
# @param kraken is the kraken file containing info on all the listed accessions

import csv
import sys
import os
#import wget
#import gzip
#import shutil

out = sys.argv[1]
kraken = sys.argv[2]
id = out[5:-23]
with open(kraken, newline="") as krakenF:
    tsv = csv.reader(krakenF, dialect=csv.excel_tab)
    firstLine = next(tsv)
    idIndex = firstLine.index("assembly_accession")
    ftpIndex = firstLine.index("ftp_path")
    for line in tsv:
        partialId = id.split("_")[0] + "_" + id.split("_")[1]
        if line[idIndex] == partialId:
            url = line[ftpIndex] + "/" + id + "_cds_from_genomic.fna.gz"
            # TODO: test for existence of wget and gzip in environment
            # Would benefit from containerization
            os.system("wget " + url + " -P ncbi/")
            os.system("gzip -d ncbi/" + id + "_cds_from_genomic.fna.gz")
            os.system("mv ncbi/" + id + "_cds_from_genomic.fna ncbi/" + id + "_cds_from_genomic.fasta")
            break