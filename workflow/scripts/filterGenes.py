# Utils for parsing fasta files
import collections
import re
import csv

# Gets the taxon id (NCBI) for the given genome id from the assembly report
# DEPRECATED: Use get_txid() to get txid from the input kraken report
# @param id is the genome id from NCBI
# @return is the taxon id
def get_txid_from_assembly(id: str) -> str:
    assembly_report = open("ncbi/" + id + "_assembly_report.txt")
    lines = assembly_report.readlines()
    for line in lines:
        if "# Taxid:" in line:
            return re.findall(r'\d+', line)[0]

# Gets the taxon id (NCBI) for the given id from the kraken report
# @param id is the genome id from NCBI
# @param kraken is the path to the kraken report
# @return is the taxon id
def get_txid(id: str, kraken: str) -> str:
    with open(kraken, newline="") as krakenF:
        tsv = csv.reader(krakenF, dialect=csv.excel_tab)
        firstLine = next(tsv)
        idIndex = firstLine.index("assembly_accession")
        txIndex = firstLine.index("taxid")
        for line in tsv:
            partialId = id.split("_")[0] + "_" + id.split("_")[1]
            if line[idIndex] == partialId:
                return line[txIndex]

# Filters sequences in a fasta file for a specific input gene
# N.B. This will only print the first occurunce of a gene in the file
# @param seqsFile is the path to the file to be parsed
# @param seq_genes is a list of genes to be filtered for (N.B. I've only tested this with list size 1)
# @return is a list containing the gene descriptor and the gene sequence
def filter_seq_genes(seqsFile: str, seq_genes: str, kraken: str) -> list:
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
                retVal = ">" + get_txid(id, kraken) + " " + " ".join(gene[0].split(" ")[1:])
                print(retVal)
                print(gene[1])
                return [retVal, gene[1]]
                break

# A wrapper for filter_seq_genes to be called from a snakemake rule
# @param seqsFile is the path to the file to be parsed
# @param output is the output file path for the snakemake rule, it is intended to take the form
#        dir/GENE__restOfFileName.fasta
# @return is undefined, instead it writes to output
def filter_seq_genes_sm(seqsFile: str, output: str, kraken: str) -> None:
    seq_genes = [output.split("/")[1].split("__")[0]]
    with open(output, "w") as out:
        vals = filter_seq_genes(seqsFile, seq_genes, kraken)
        for val in vals:
            out.write(val + "\n")