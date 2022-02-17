# Downloads genome from NCBI for a given genome id
import csv
import os

# Downloads the latest complete genome matching the input genome id (which has a 1-to-1 correspondence with a species-level taxon id)
# @param id is the genome id
# @param ext is the file extension (specifying cds or rna)
# @return is the exit code from the last os.system call
def download_genome(genomeId: str, ext: str) -> str:
    with open("run_assembly.txt") as assembly:
        tsv = csv.reader(assembly, dialect=csv.excel_tab)
        firstLine = next(tsv)
        accIndex = firstLine.index("assembly_accession")
        ftpIndex = firstLine.index("ftp_path")

        for line in tsv:
            if line[accIndex] == (genomeId.split("_")[0] + "_" + genomeId.split("_")[1]):
                url = line[ftpIndex] + "/" + genomeId + ext
                #url = "rsync" + line[ftpIndex][5:] + "/" + line[accIndex] + ext
                os.system("wget " + url + " -P ncbi/")
                #os.system("rsync --copy-links --times --verbose " + url + " ncbi/") # Recommended method but weird
                os.system("gzip -d ncbi/" + genomeId + ext)
                return os.system("mv ncbi/" + genomeId + ext[:-3] + " ncbi/" + genomeId + ext[:-6] + "fasta")

# Wrapper for download_genes from snakemake
# @param output is the output file containing the genome id
# @return is undefined
def download_genes_sm(output: str) -> None:
    outputs = output.split(" ")
    for out in outputs:
        genomeId = out[5:-23]
        ext = out[-23:-5] + "fna.gz"
        download_genome(genomeId, ext)














