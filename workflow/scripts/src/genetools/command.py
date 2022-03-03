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
    prepare_run_assembly(args.taxon_list, args.output_dir, logF)
    downloaded_genome_ids, failed_genome_ids = download_genomes(args.output_dir, logF)

    if failed_genome_ids:
        logF.write("Failed to download genome(s) for:\n")
        for failure in failed_genome_ids:
            logF.write(failure + "\n")
    
    return downloaded_genome_ids

def _filter_genes(args, downloaded_genome_ids, logF):
    gene_counter = extract_genes(args.filter_genes.name, downloaded_genome_ids, args.output_dir)

    logF.write("Gene extraction successes (total=" + str(len(downloaded_genome_ids)) + "):\n")
    for gene in gene_counter.most_common():
        logF.write(gene[0] + ": " + str(gene[1]) + "\n")
    logF.write("Ignoring genes with fewer than 4 sequences...\n")

    return gene_counter

def _merge_genes(args, gene_counter):
    merge_genes(gene_counter, args.output_dir)

def _write_config(args, cfg_fp):
    print("Writing " + cfg_fp + " for this run")
    with open(cfg_fp, "w") as config:
        config.write("# Config file for tree building pipeline\n")
        config.write("# Genes to build trees from\n")
        config.write("GENES: [")
        for gene in args.filter_genes.readlines():
            if gene[0] != "#":
                config.write("\"" + gene.strip() + "\", ")
        config.write("]")

def dir_path(string):
    if string == "":
        return ""
    elif os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("taxon_list", help="Filepath to list of taxa to download genomes for (tsv with taxon id in first column)")
    p.add_argument("-f", "--filter_genes", type=argparse.FileType("r"),
        help=(
            "Filepath to list of genes to filter into separate files (one gene per line)"))
    p.add_argument("-m", "--merge_genes", action="store_true",
        help=(
            "Filepath of metadata. Example is example_metadata.txt"))
    p.add_argument("--remove_temp", action="store_true",
        help=(
            "Removes any files that aren't produced in the final step of the program (i.e. if run with --merge-genes, this would remove the downloaded sequences and individual gene files, keeping only the merged gene files)"))
    p.add_argument("--write_config", action="store_true",
        help=(
            "Writes a config file for snakemake based on the genes provided in the --filter_genes arg"))
    p.add_argument("-o", "--output_dir", default="", type=dir_path,
        help=(
            "Location to write outputs to, can be an absolute path or a path relative to the library base"))

    args = p.parse_args(argv)

    if args.merge_genes and not args.filter_genes:
        p.error("Must include --filter_genes arg along with --merge_genes")
    
    if args.write_config and not args.filter_genes:
        p.error("Must include --filter_genes arg along with --write_config")
    
    None if os.path.isdir(os.path.join(args.output_dir, "output/")) else os.mkdir(os.path.join(args.output_dir, "output/"))
    None if os.path.isdir(os.path.join(args.output_dir, "output/ncbi/")) else os.mkdir(os.path.join(args.output_dir, "output/ncbi"))
    None if os.path.isdir(os.path.join(args.output_dir, "logs/")) else os.mkdir(os.path.join(args.output_dir, "logs"))
    logF = open(os.path.join(args.output_dir, "logs/main.log"), "w")
    print("Writing logs to: " + logF.name)

    logF.write("Downloading genomes...\n")
    ids = _download_genes(args, logF)

    gene_counter = Counter()
    if args.filter_genes:
        None if os.path.isdir(os.path.join(args.output_dir, "output/sequences/")) else os.mkdir(os.path.join(args.output_dir, "output/sequences"))
        None if os.path.isdir(os.path.join(args.output_dir, "data/")) else os.mkdir(os.path.join(args.output_dir, "data"))
        logF.write("Extracting genes...\n")
        gene_counter = _filter_genes(args, ids, logF)
        if args.remove_temp:
            shutil.rmtree(os.path.join(args.output_dir, "output/ncbi/"))
    
    if args.merge_genes:
        None if os.path.isdir(os.path.join(args.output_dir, "output/merged-sequences/")) else os.mkdir(os.path.join(args.output_dir, "output/merged-sequences"))
        logF.write("Merging gene files...\n")
        _merge_genes(args, gene_counter)
        if args.remove_temp:
            shutil.rmtree(os.path.join(args.output_dir, "output/sequences/"))
    
    if args.write_config:
        cfg_fp = os.path.join(args.output_dir, "output/config.yml")
        logF.write("Writing config to " + cfg_fp)
        _write_config(args, cfg_fp)

    logF.close()

    