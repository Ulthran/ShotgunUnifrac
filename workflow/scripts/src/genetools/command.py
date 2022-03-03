import argparse
import csv
import os
import shutil
import sys
import tqdm
from collections import Counter
from typing import Tuple

from .downloadGenes import prepare_run_assembly, download_genomes
from .filterGenes import extract_genes
from .mergeGenes import merge_genes

def _download_genes(args, logF):
    prepare_run_assembly(args.taxon_list, logF)
    downloaded_genome_ids, failed_genome_ids = download_genomes(logF)

    if failed_genome_ids:
        logF.write("Failed to download genome(s) for:\n")
        for failure in failed_genome_ids:
            logF.write(failure + "\n")
    
    return downloaded_genome_ids

def _filter_genes(args, downloaded_genome_ids, logF):
    gene_counter = extract_genes(args.filter_genes.name, downloaded_genome_ids)

    logF.write("Gene extraction successes (total=" + str(len(downloaded_genome_ids)) + "):\n")
    for gene in gene_counter.most_common():
        logF.write(gene[0] + ": " + str(gene[1]) + "\n")
    logF.write("Ignoring genes with fewer than 4 sequences...\n")

    return gene_counter

def _merge_genes(gene_counter):
    merge_genes(gene_counter)

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("taxon_list", help="Filepath to list of taxa to download genomes for (tsv with taxon id in first column)")
    p.add_argument("--filter_genes", type=argparse.FileType("r"),
        help=(
            "Filepath to list of genes to filter into separate files (one gene per line)"))
    p.add_argument("--merge_genes", action="store_true",
        help=(
            "Filepath of metadata. Example is example_metadata.txt"))
    p.add_argument("--remove_temp", action="store_true",
        help=(
            "Removes any files that aren't produced in the final step of the program (i.e. if run with --merge-genes, this would remove the downloaded sequences and individual gene files, keeping only the merged gene files)"))
    
    args = p.parse_args(argv)

    if args.merge_genes and not args.filter_genes:
        p.error("Must include --filter_genes arg along with --merge_genes")
    
    None if os.path.isdir("output/") else os.mkdir("output")
    None if os.path.isdir("output/ncbi/") else os.mkdir("output/ncbi")
    None if os.path.isdir("logs/") else os.mkdir("logs")
    logF = open("logs/main.log", "w")
    print("Writing logs to: " + logF.name)

    logF.write("Downloading genomes...\n")
    ids = _download_genes(args, logF)

    gene_counter = Counter()
    if args.filter_genes:
        None if os.path.isdir("output/sequences/") else os.mkdir("output/sequences")
        logF.write("Extracting genes...\n")
        gene_counter = _filter_genes(args, ids, logF)
        if args.remove_temp:
            shutil.rmtree("output/ncbi/")
    
    if args.merge_genes:
        None if os.path.isdir("output/merged-sequences/") else os.mkdir("output/merged-sequences")
        logF.write("Merging gene files...\n")
        _merge_genes(gene_counter)
        if args.remove_temp:
            shutil.rmtree("output/sequences/")



    