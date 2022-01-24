

import collections
import re

def get_seq_gene(desc):
    return desc.split()[1]

def filter_seq_genes(seqsFile, seq_genes):
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
                print(gene[0])
                print(gene[1])
                return gene[0], gene[1]

def filter_seq_genes_sm(seqsFile, output):
    seq_genes = [output.split("/")[1].split("__")[0]]
    return filter_seq_genes(seqsFile, seq_genes)