# Merge individual gene files
import os
from collections import Counter

# Merges extracted genes together into single files
# @param gene_counter is the Counter containing counts of all extraced genes
# @param outputDir is the location for merged gene files
# @return is a list of genes with more than 3 counts (enough to contruct a tree)
def merge_genes(gene_counter: Counter, outputDir: str) -> list:
    seq_dir = os.path.join(outputDir, "output/sequences/")
    mer_seq_dir = os.path.join(outputDir, "output/merged-sequences/")
    cogs = set()

    for f in os.listdir(os.fsencode(seq_dir)):
        fn = os.fsdecode(f)
        cogs.add(fn.split("__")[0])

    
    for cog in cogs:
        with open(mer_seq_dir + cog + ".fasta", "w") as mergedF:
            for f in os.listdir(os.fsencode(seq_dir)):
                fn = os.fsdecode(f)
                if cog in fn:
                    with open(seq_dir + fn) as seqF:
                        for l in seqF.readlines():
                            mergedF.write(l)
    return cogs