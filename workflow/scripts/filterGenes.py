# Utils for parsing fasta files
import collections
import re
import csv
from typing import Union
import sys

# Gets the taxon id (NCBI) for the given id from the kraken report
# @param id is the genome id from NCBI
# @param run_path is the path to run_assembly.txt
# @return is the taxon id
def get_txid(id: str, run_path: str = "") -> str:
    with open(run_path + "run_assembly.txt") as runAssembly:
        tsv = csv.reader(runAssembly, dialect=csv.excel_tab)
        firstLine = next(tsv)
        idIndex = firstLine.index("assembly_accession")
        txIndex = firstLine.index("species_taxid")
        partialId = id.split("_")[0] + "_" + id.split("_")[1]
        for line in tsv:
            if line[idIndex] == partialId:
                return line[txIndex]

# Filter through the given fasta file and find targetGene if it's there
# @param seqFile is the path to the fasta file
# @param genomeId is the id of the desired genome
# @param targetGene is the gene to be found
# @return is either a list containing the description/gene pair or None if it couldn't be found
def filter_seq_file(seqFile: str, genomeId: str, targetGene: str) -> Union[list, None]:
    run_assembly_path = "" if len(seqFile.split("ncbi")) == 1 else seqFile.split("ncbi")[0]
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
    for gene in seqList:
        # Search for [gene=GENENAME] pattern
        seq_gene = gene[0].split()[1] # If there is a gene tag, it should be the fist item after the id
        match = re.search("\[gene=.*\]", seq_gene)
        if match:
            gene_name = seq_gene.strip('[gene=').strip(']')
            if gene_name.upper() == targetGene.upper():
                retVal = ">" + get_txid(genomeId, run_assembly_path) + " " + " ".join(gene[0].split(" ")[1:])
                #print(retVal)
                #print(gene[1])
                return [retVal, gene[1]]
        
        # Search for [protein=*GENENAME*] pattern
        protIndex = gene[0].find("[protein=")
        endProtIndex = gene[0].find("]", protIndex)
        proteinStr = gene[0][protIndex:endProtIndex]
        if targetGene.upper() in proteinStr.upper():
            retVal = ">" + get_txid(genomeId, run_assembly_path) + " " + " ".join(gene[0].split(" ")[1:])
            #print(retVal)
            #print(gene[1])
            return [retVal, gene[1]]

        # Search for [product=*GENENAME*] pattern
        prodIndex = gene[0].find("[product=")
        endProdIndex = gene[0].find("]", prodIndex)
        productStr = gene[0][prodIndex:endProdIndex]
        if targetGene.upper() in productStr.upper():
            retVal = ">" + get_txid(genomeId, run_assembly_path) + " " + " ".join(gene[0].split(" ")[1:])
            #print(retVal)
            #print(gene[1])
            return [retVal, gene[1]]
    
    return None


# Filters sequences in a fasta file for a specific input gene
# N.B. This will only find the first occurunce of a gene in the file
# @param genomeId is the genome id for the fasta file to be parsed
# @param gene is the gene to be found
# @param ncbi_dir is the directory containing downloaded genomes
# @return is a list containing the gene descriptor and the gene sequence
def filter_seq_genes(genomeId: str, gene: str, ncbi_dir: str) -> list:
    cds_path = ncbi_dir + genomeId + "_cds_from_genomic.fasta"
    rna_path = ncbi_dir + genomeId + "_rna_from_genomic.fasta"

    vals = filter_seq_file(cds_path, genomeId, gene)
    if not vals:
        vals = filter_seq_file(rna_path, genomeId, gene)
    if not vals:
        return []
    return vals

# A wrapper for filter_seq_genes to be called from a snakemake rule
# @param seqsFile is the path to the file to be parsed
# @param output is the output file path for the snakemake rule, it is intended to take the form
#        dir/GENE__restOfFileName.fasta
# @param genomeId is the genome id that can be passed if it's there when calling to increase the robustness of this function
# @return is None, instead it writes to output
def filter_seq_genes_sm(genomeId: str, gene: str, ncbi_dir: str, seq_dir: str) -> None:
    with open(seq_dir + gene + "__" + genomeId + ".fasta", "w") as out:
        try:
            vals = filter_seq_genes(genomeId, gene, ncbi_dir)
            for val in vals:
                out.write(val + "\n")
        except TypeError:
            None