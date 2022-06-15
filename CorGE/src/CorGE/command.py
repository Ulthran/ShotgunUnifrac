import argparse
import logging
from pathlib import Path
from .collect import collect_genomes
from .extract import extract_genes

def collect_genomes(args: argparse.Namespace):
    logging.debug(args)
    ncbi_species = open(args.ncbi_species) if args.ncbi_species else None
    ncbi_accessions = open(args.ncbi_accessions) if args.ncbi_accessions else None
    collect_genomes(args.output_dir, ncbi_species, ncbi_accessions, args.local, args.outgroup)
    ncbi_species.close()
    ncbi_accessions.close()

def extract_genes(args: argparse.Namespace):
    logging.debug(args)
    extract_genes(args.genomes, args.output)

def dir_path(path):
    if Path.is_dir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

def main(argv=None):
    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(help='Subcommands')

    collect_genomes_subparser = subparsers.add_parser(
        "collect_genomes",
        help="Collect nucleotide- and protein-encoded genomes of interest")
    extract_genes_subparser = subparsers.add_parser(
        "extract_genes",
        help="Extract SCCGs from all collected genomes and curate data for tree building")

    collect_genomes_subparser.add_argument("output_dir",
        type=dir_path,
        help="Directory to collect genomes in")
    collect_genomes_subparser.add_argument("--ncbi_species",
        type=argparse.FileType("r"),
        help="File listing species level taxon ids to be collected from NCBI")
    collect_genomes_subparser.add_argument("--ncbi_accessions",
        type=argparse.FileType("r"),
        help="File listing genome accessions to be collected from NCBI")
    collect_genomes_subparser.add_argument("--local",
        type=dir_path,
        default="",
        help="Directory containing nucleotide- and protein-encoded pairs of genome files. Any unpaired files will be ignored")
    collect_genomes_subparser.add_argument("--outgroup",
        type=str,
        default="2173",
        help="Specify the outgroup for tree rooting. Integers will be parsed as species level taxon ids and retrieved from NCBI. Otherwise will search for a matching nucleotide-encoded file in ouput_dir or local (Default: 2173, enter None to not use outgroup rooting)")
    collect_genomes_subparser.set_defaults(func=collect_genomes)

    extract_genes_subparser.add_argument("genomes",
        type=dir_path,
        help="Directory with collected genomes (curated with collect_genomes)")
    extract_genes_subparser.add_argument("-o", "--output",
        type=dir_path,
        default="./",
        help="Directory to write output to (Default: ./)")
    extract_genes_subparser.set_defaults(func=extract_genes)

    args = main_parser.parse_args(argv)

    args.func(args)

