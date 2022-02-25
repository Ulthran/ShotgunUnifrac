# Downloads genome from NCBI for a given genome id
import csv
import os
import tqdm
from typing import Tuple
from collections import Counter
import filterGenes

# Downloads genomes from NCBI for all the accessions in run_assembly.txt
# @param prefix is a prefix to the path (for testing)
# @return is a tuple containing the ids of successfully downloaded genomes and those of failed genomes
def download_genomes(prefix: str = "") -> Tuple[list, list]:
    downloaded_genome_ids = []
    failed_genome_ids = []
    ncbi_dir = "ncbi/" if prefix == "" else prefix + "/ncbi/"
    run_assembly_path = "run_assembly.txt" if prefix == "" else prefix + "/run_assembly.txt"
    with open (run_assembly_path) as run_assembly:
        run_assembly_reader = csv.reader(run_assembly, dialect=csv.excel_tab)
        firstLine = next(run_assembly_reader)
        ftpIndex = firstLine.index("ftp_path")

        for line in run_assembly_reader:
            genomeId = line[ftpIndex].split("/")[-1]
            for ext in ["_cds_from_genomic.fna.gz", "_rna_from_genomic.fna.gz"]:
                if not os.path.isfile(ncbi_dir + genomeId + ext[:-6] + "fasta"):
                    url = line[ftpIndex] + "/" + genomeId + ext
                    os.system("wget " + url + " -P " + ncbi_dir)
                    os.system("gzip -d " + ncbi_dir + genomeId + ext)
                    os.system("mv " + ncbi_dir + genomeId + ext[:-3] + " " + ncbi_dir + genomeId + ext[:-6] + "fasta")
            if os.path.isfile(ncbi_dir + genomeId + "_cds_from_genomic.fasta") and os.path.isfile(ncbi_dir + genomeId + "_rna_from_genomic.fasta"):
                downloaded_genome_ids.append(genomeId)
            else:
                failed_genome_ids.append(genomeId)
    return downloaded_genome_ids, failed_genome_ids

# Extracts each gene in geneFile from each genome in downloaded_genome_ids
# @param geneFile is the path to the file listing genes
# @param downloaded_genome_ids is the list of successfully downloaded genome ids
# @param prefix is a prefix to the path (for testing)
# @return is a Counter containing how many of each gene was successfully extracted
def extract_genes(geneFile: str, downloaded_genome_ids: list, prefix: str = "") -> Counter:
    gene_counter = Counter()
    with open(geneFile) as gene_file: # Initialize Counter so that it picks up genes that don't exist in any genomes
        gene_file_reader = csv.reader(gene_file)
        for line in gene_file_reader:
            if line[0][0] != "#":
                gene_counter[line[0]] = 0

    ncbi_dir = "ncbi/" if prefix == "" else prefix + "/ncbi/"
    seq_dir = "sequences/" if prefix == "" else prefix + "/sequences/"
    print("Extracting genes from genome files...\n")
    print(str(downloaded_genome_ids))
    with tqdm.tqdm(total=len(downloaded_genome_ids) * len(gene_counter.most_common())) as pbar: # Add progress bar
        for genomeId in downloaded_genome_ids:
            with open(geneFile) as gene_file:
                gene_file_reader = csv.reader(gene_file)
                for line in gene_file_reader:
                    if line[0][0] != "#":
                        filterGenes.filter_seq_genes_sm(
                                genomeId,
                                line[0],
                                ncbi_dir,
                                seq_dir
                            )
                        pbar.update(1)
                        if os.path.getsize(seq_dir + line[0] + "__" + genomeId + ".fasta") != 0:
                            gene_counter[line[0]] += 1
    return gene_counter

# Merges extracted genes together into single files
# @param gene_counter is the Counter containing counts of all extraced genes
# @param prefix is a prefix to the path (for testing)
# @return is a list of genes with more than 3 counts (enough to contruct a tree)
def merge_genes(gene_counter: Counter, prefix: str = "") -> list:
    gene_list = []
    for gene in gene_counter.most_common():
        if gene[1] > 3: # Minimum number of taxa to infer evolutionary history
            gene_list.append(gene[0])

    seq_dir = "sequences/" if prefix == "" else prefix + "/sequences/"
    mer_seq_dir = "merged-sequences/" if prefix == "" else prefix + "/merged-sequences/"
    for gene in gene_list:
        with open(mer_seq_dir + gene + ".fasta", "w") as mergedF:
            directory = os.fsencode(seq_dir)
            for file in os.listdir(directory):
                filename = os.fsdecode(file)
                if gene in filename: 
                    with open(seq_dir + filename) as seqF:
                        for line in seqF.readlines():
                            mergedF.write(line)
    return gene_list

# Wrapper for download_genes from snakemake
# @param output is the output file containing the genome id
# @return is undefined
def download_genes_sm(output: str) -> None:
    outputs = output.split(" ")
    for out in outputs:
        genomeId = out[5:-23]
        ext = out[-23:-5] + "fna.gz"
        download_genome(genomeId, ext)














