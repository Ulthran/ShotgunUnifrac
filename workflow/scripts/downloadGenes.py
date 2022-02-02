# Downloads genome from NCBI for a given genome id
# @param out is the output file continaing the accession of the desired genome for download
#        taking the form ncbi/{ID}_cds_from_genomic.fasta
# @param kraken is the kraken file containing info on all the listed accessions

import csv
import sys
import os

# Downloads genome from NCBI
# @param id is the genome id
# @param kraken is the path to the kraken report containing the FTP path
# @return is undefined, os.system calls do the work
def download_genes(id: str, kraken: str) -> None:
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

# Wrapper for download_genes from snakemake
# @param output is the output file containing the genome id
# @param kraken is the kraken report file path
# @return is undefined
def download_genes_sm(output: str, kraken: str) -> None:
    id = output[5:-23]
    download_genes(id, kraken)