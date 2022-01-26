import collections
import re

# Gets the taxon id (NCBI) for the given genome id from the assembly report
# @param id is the genome id from NCBI
# @return is the taxon id
def get_txid(id):
    assembly_report = open("ncbi/" + id + "_assembly_report.txt")
    lines = assembly_report.readlines()
    for line in lines:
        if "# Taxid:" in line:
            return re.findall(r'\d+', line)[0]

# Filters sequences in a fasta file for a specific input gene
# N.B. This will only print the first occurunce of a gene in the file
# @param seqsFile is the path to the file to be parsed
# @param seq_genes is a list of genes to be filtered for (N.B. I've only tested this with list size 1)
# @return is undefined, this is intended to be run from the shell and then stdout is captured
def filter_seq_genes(seqsFile, seq_genes):
    id = seqsFile.split("/")[1].split("_cds")[0]
    fasta = open(seqsFile, "r")
    seqs = fasta.readlines()
    seqList = []
    seqObj = []

    for obj in seqs:
        if(obj[0] == ">"):
            seqList.append(seqObj)
            seqObj = ["", ""]
            seqObj[0] = obj.rstrip()
        else:
            seqObj[1] = seqObj[1] + obj.rstrip()
    seqList.pop(0) # Get rid of blank entry at the start

    for gene in seqList:
        seq_gene = gene[0].split()[1]
        match = re.search("\[gene=.*\]", seq_gene)
        if match:
            gene_name = seq_gene.strip('[gene=').strip(']')
            if gene_name in seq_genes:
                retVal = ">" + get_txid(id) + " " + " ".join(gene[0].split(" ")[1:])
                print(retVal)
                print(gene[1])
                break

# A wrapper for filter_seq_genes to be called from a snakemake rule
# @param seqsFile is the path to the file to be parsed
# @param output is the output file path for the snakemake rule, it is intended to take the form
#        dir/GENE__restOfFileName.fasta
# @return is undefined, see filter_seq_genes
def filter_seq_genes_sm(seqsFile, output):
    seq_genes = [output.split("/")[1].split("__")[0]]
    filter_seq_genes(seqsFile, seq_genes)