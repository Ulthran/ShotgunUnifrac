import collections
import re

def get_seq_gene(desc):
    return desc.split()[1]

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
        seq_gene = get_seq_gene(gene[0])
        match = re.search("\[gene=.*\]", seq_gene)
        if match:
            gene_name = seq_gene.strip('[gene=').strip(']')
            if gene_name in seq_genes:
                retVal = ">" + id + " " + " ".join(gene[0].split(" ")[1:])
                print(retVal)
                print(gene[1])

def filter_seq_genes_sm(seqsFile, output):
    seq_genes = [output.split("/")[1].split("__")[0]]
    filter_seq_genes(seqsFile, seq_genes)