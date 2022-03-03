# Merge individual gene files
import os
from collections import Counter

# Merges extracted genes together into single files
# @param gene_counter is the Counter containing counts of all extraced genes
# @param outputDir is the location for merged gene files
# @return is a list of genes with more than 3 counts (enough to contruct a tree)
def merge_genes(gene_counter: Counter, outputDir: str) -> list:
    gene_list = []
    for gene in gene_counter.most_common():
        if gene[1] > 3: # Minimum number of taxa to infer evolutionary history
            gene_list.append(gene[0])

    seq_dir = os.path.join(outputDir, "output/sequences/")
    mer_seq_dir = os.path.join(outputDir, "output/merged-sequences/")
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