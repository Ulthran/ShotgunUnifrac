# Utils for parsing fasta files
import collections
from io import TextIOWrapper
import re
import csv
from typing import Union
import sys
import os
from collections import Counter
import tqdm

# Gets the taxon id (NCBI) for the given id from the kraken report
# @param id is the genome id from NCBI
# @param run_path is the path to run_assembly.txt
# @return is the taxon id
def get_txid(id: str, run_path: str = "") -> str:
    with open(os.path.join(run_path, "data/run_assembly.txt")) as runAssembly:
        tsv = csv.reader(runAssembly, dialect=csv.excel_tab)
        firstLine = next(tsv)
        idIndex = firstLine.index("assembly_accession")
        txIndex = firstLine.index("species_taxid")
        partialId = id.split("_")[0] + "_" + id.split("_")[1]
        for line in tsv:
            if line[idIndex] == partialId:
                return line[txIndex]

# Finds the given pattern in the gene description
# @param pattern is the pattern to search for (i.e. "[gene=")
# @param description is the gene's description to pattern match on
# @param endPattern is the end of the pattern to search for (i.e. "]")
# @return is the uppercase search result
def find_gene(pattern: str, description: str, endPattern: str = "]") -> str:
    # Search for [pattern=*GENENAME*]
    geneIndex = description.find(pattern)
    endGeneIndex = description.find(endPattern, geneIndex)
    geneStr = description[geneIndex:endGeneIndex]
    return geneStr.upper()

# Finds genes not found by annotation using vsearch
# @param refDB is the path to the reference DB fasta file
# @param queryDB is the path to the query fasta
# @param out is the path to the output file, should be in the proper output location for gene sequences (i.e. output/sequences/)
# @param txid is the taxon id of the given refDB fasta
# @return is True if the gene was found, False otherwise
def vsearch_gene(refDB: str, queryDB: str, out: str, txid: str) -> bool:
    # Ex: vsearch --usearch_global output/ncbi/GCF_000008865.2_ASM886v2_cds_from_genomic.fasta --db secG.fasta --fastapairs out.out --id 0.9
    os.system(f"vsearch --usearch_global {refDB} --db {queryDB} --fastapairs {out} --id 0.9")
    with open(out) as f:
        data = f.readlines()
    if data != []:
        with open(out, "w") as f:
            f.write(f"> {txid}\n")
            for line in data[1:]:
                if line[0] != ">":
                    f.write(line)
                else:
                    break
        return True
    else:
        return False

# Filter through the given fasta file and find targetGene if it's there
# @param seqFile is the path to the fasta file
# @param genomeId is the id of the genome in seqFile (needed for getting the taxon id)
# @param targetGene is the gene to be found
# @param rename is a boolean telling the function whether or not to reannotate gene descriptions
# @return is either a list containing the description/gene pair or None if it was written with vsearch or couldn't be found
def filter_seq_file(seqFile: str, genomeId: str, targetGene: str, rename: bool = True) -> Union[list, None]:
    run_assembly_path = "" if len(seqFile.split("output/")) == 1 else seqFile.split("output/")[0]
    seqList = []
    seqObj = []
    # Fill seqList with all the description/sequence pairs from the file
    with open(seqFile) as seqF:
        # Loop through file to create list of description/sequence pairs of each gene
        for obj in seqF.readlines():
            if(obj[0] == ">"):
                seqList.append(seqObj)
                seqObj = ["", ""]
                seqObj[0] = obj.rstrip()
            else:
                seqObj[1] = seqObj[1] + obj.rstrip()
        seqList.pop(0) # Get rid of blank entry at the start
        seqList.append(seqObj)
    
    # Find the desired gene and return the description/sequence pair
    # Searching [gene=*], [protein=*], and [product=*]
    for gene in seqList:
        if targetGene.upper() in find_gene("[gene=", gene[0]):
            retVal = ">" + get_txid(genomeId, run_assembly_path) if rename else gene[0]
            return [retVal, gene[1]]
        
        if targetGene.upper() in find_gene("[protein=", gene[0]):
            retVal = ">" + get_txid(genomeId, run_assembly_path) if rename else gene[0]
            return [retVal, gene[1]]

        if targetGene.upper() in find_gene("[product=", gene[0]):
            retVal = ">" + get_txid(genomeId, run_assembly_path) if rename else gene[0]
            return [retVal, gene[1]]
    
    # Search with vsearch if annotation pattern matching didn't find anything
    target_fp = os.path.join(run_assembly_path, f"output/ref-genes/{targetGene}.fasta")
    output_fp = os.path.join(run_assembly_path, f"output/sequences/{targetGene}__{genomeId}.fasta")
    if vsearch_gene(seqFile, target_fp, output_fp, get_txid(genomeId, run_assembly_path)):
        print(f"Found {targetGene} in {seqFile} using vsearch")
        return None
    else:
        print(f"Didn't find {targetGene} in {seqFile} using vsearch")
    
    return None


# Filters sequences in a fasta file for a specific input gene
# N.B. This will only find the first occurunce of a gene in the file
# @param genomeId is the genome id for the fasta file to be parsed
# @param gene is the gene to be found
# @param ncbi_dir is the directory containing downloaded genomes
# @param rename is a boolean telling the function whether or not to reannotate gene descriptions
# @return is a list containing the gene descriptor and the gene sequence
def filter_seq_genes(genomeId: str, gene: str, ncbi_dir: str, rename: bool = True) -> list:
    directory = os.fsencode(ncbi_dir)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if (genomeId in filename and (filename.split('.')[-1] == "fasta" or filename.split('.')[-1] == "fna")
            and (("tRNA" in gene and "rna" in filename) or ("tRNA" not in gene and "cds" in filename))):
            vals = filter_seq_file(os.path.join(ncbi_dir, filename), genomeId, gene, rename)
            if vals:
                return vals
    return []

# Extracts each gene in geneFile from each genome in downloaded_genome_ids
# @param geneFile is the path to the file listing genes
# @param downloaded_genome_ids is the list of successfully downloaded genome ids
# @param outputDir is the location to extract genes to
# @param inputDir is an alternative to the default ncbi directory to find genomes in
# @return is a Counter containing how many of each gene was successfully extracted
def extract_genes(geneFile: str, downloaded_genome_ids: list, outputDir: str, inputDir: str = "") -> Counter:
    gene_counter = Counter()
    with open(geneFile) as gene_file: # Initialize Counter so that it picks up genes that don't exist in any genomes
        gene_file_reader = csv.reader(gene_file)
        for line in gene_file_reader:
            if line[0][0] != "#":
                gene_counter[line[0]] = 0

    ncbi_dir = os.path.join(outputDir, "output/ncbi/")
    if inputDir != "": # Handle alternative input, add each file in dir to downloaded_genome_ids
        directory = os.fsencode(inputDir)
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            downloaded_genome_ids.append('.'.join(filename.split(".")[:-1]))
        ncbi_dir = inputDir
    seq_dir = os.path.join(outputDir, "output/sequences/")
    print("Extracting genes from genome files...\n")
    print(str(downloaded_genome_ids))
    with tqdm.tqdm(total=len(downloaded_genome_ids) * len(gene_counter.most_common())) as pbar: # Add progress bar
        for genomeId in downloaded_genome_ids:
            with open(geneFile) as gene_file:
                gene_file_reader = csv.reader(gene_file)
                for line in gene_file_reader:
                    if line[0][0] != "#":
                        with open(seq_dir + line[0] + "__" + genomeId + ".fasta", "w") as out:
                            try:
                                vals = []
                                if inputDir != "":
                                    vals = filter_seq_genes(genomeId, line[0], ncbi_dir, rename=False)
                                else:
                                    vals = filter_seq_genes(genomeId, line[0], ncbi_dir, rename=True)
                                for val in vals:
                                    out.write(val + "\n")
                            except TypeError:
                                None
                        pbar.update(1)
                        if os.path.getsize(seq_dir + line[0] + "__" + genomeId + ".fasta") != 0:
                            gene_counter[line[0]] += 1
    return gene_counter
