import argparse
import logging
import os
import sys
from enum import Enum
from pathlib import Path
from .collect import collect_genomes
from .extract import extract_genes


class FileType(Enum):
    prot = "prot"
    nucl = "nucl"

    def __str__(self):
        return self.value


class NameType(Enum):
    acc = "acc"
    txid = "txid"
    strain = "strain"
    species = "species"

    def __str__(self):
        return self.value


def _collect_genomes(args: argparse.Namespace):
    if args.all and (args.ncbi_species or args.ncbi_accessions or args.local_fp):
        logging.error(
            "Can't use --all with --ncbi_species, --ncbi_accessions, or --local"
        )
        sys.exit(1)

    logging.getLogger().setLevel(args.log_level)
    collect_genomes(vars(args))


def _extract_genes(args: argparse.Namespace):
    logging.getLogger().setLevel(args.log_level)
    extract_genes(vars(args))


def dir_path(dir: str):
    if Path.is_dir(Path(dir)):
        return dir
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{dir} is not a valid path")


def main(argv=None):
    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(help="Subcommands")

    collect_genomes_subparser = subparsers.add_parser(
        "collect_genomes",
        help="Collect nucleotide- and protein-encoded genomes of interest",
    )
    extract_genes_subparser = subparsers.add_parser(
        "extract_genes",
        help="Extract SCCGs from all collected genomes and curate data for tree building",
    )

    collect_genomes_subparser.add_argument(
        "--output_fp",
        type=dir_path,
        help="Directory to collect genomes in (Default: ./output/)",
    )
    collect_genomes_subparser.add_argument(
        "--all",
        action="store_true",
        default=False,
        help="Collect one representative genome from each species listed in NCBI's RefSeq database. Don't use this with --ncbi_species, --ncbi_accessions, or --local",
    )
    collect_genomes_subparser.add_argument(
        "--ncbi_species",
        type=argparse.FileType("r"),
        help="File listing species level taxon ids to be collected from NCBI",
    )
    collect_genomes_subparser.add_argument(
        "--ncbi_accessions",
        type=argparse.FileType("r"),
        help="File listing genome accessions to be collected from NCBI",
    )
    collect_genomes_subparser.add_argument(
        "--local_fp",
        type=dir_path,
        help="Directory containing nucleotide- and protein-encoded pairs of genome files. Any unpaired files will be ignored",
    )
    collect_genomes_subparser.add_argument(
        "--outgroup",
        type=str,
        help="Specify the outgroup for tree rooting. Integers will be parsed as species level taxon ids and retrieved from NCBI. Otherwise will search for a matching nucleotide-encoded file in ouput_dir or local (Default: 2173, enter None to not use outgroup rooting)",
    )
    collect_genomes_subparser.add_argument(
        "-n",
        action="store_true",
        default=False,
        help="Dry run, show what would be gathered but don't do it",
    )
    collect_genomes_subparser.add_argument(
        "--log_level",
        type=int,
        default=20,
        help="Sets the log level, default is info, 10 for debug (Default: 20)",
    )
    collect_genomes_subparser.set_defaults(func=_collect_genomes)

    extract_genes_subparser.add_argument(
        "--genomes",
        type=dir_path,
        help="Directory with collected genomes (curated with collect_genomes) (Default: ./output/genomes/)",
    )
    extract_genes_subparser.add_argument(
        "-o",
        "--output",
        type=dir_path,
        help="Directory to write output to (Default: ./output/)",
    )
    extract_genes_subparser.add_argument(
        "-t",
        "--file_type",
        type=FileType,
        choices=list(FileType),
        help="Output in merged-sequences can be nucleotide- or protein-encoded (Default: prot)",
    )
    extract_genes_subparser.add_argument(
        "-n",
        "--name_type",
        type=NameType,
        choices=list(NameType),
        help="Names to show on final tree (Default: txid)",
    )
    extract_genes_subparser.add_argument(
        "--outgroup",
        type=str,
        help="Outgroup to use for tree rooting, name must correspond with files in the genomes dir (Default: 2173)",
    )
    extract_genes_subparser.add_argument(
        "--log_level",
        type=int,
        default=20,
        help="Sets the log level, default is info, 10 for debug (Default: 20)",
    )
    extract_genes_subparser.set_defaults(func=_extract_genes)

    args = main_parser.parse_args(argv)

    # print usage if user only enters the program name
    if len(sys.argv) < 2:
        main_parser.print_usage()
        sys.exit(1)

    logging.basicConfig()
    args.func(args)

    #return main_parser
