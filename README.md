# ShotgunUnifrac

<!--Begin status badges-->
[![CI](https://github.com/Ulthran/ShotgunUnifrac/actions/workflows/main.yml/badge.svg)](https://github.com/Ulthran/ShotgunUnifrac/actions/workflows/main.yml)
[![codecov](https://codecov.io/gh/Ulthran/ShotgunUnifrac/branch/master/graph/badge.svg?token=N9KSWRS4XG)](https://codecov.io/gh/Ulthran/ShotgunUnifrac)
[![Documentation Status](https://readthedocs.org/projects/shotgununifrac/badge/?version=latest)](https://shotgununifrac.readthedocs.io/en/latest/?badge=latest)
<!--End status badges-->

A dual use program for downloading and extracting genes from NCBI and for creating phylogenetic trees for many marker genes and merging the results into one

## Install

    git clone git@github.com:Ulthran/ShotgunUnifrac.git

To install the `CorGE` library for downloading, extracting, and merging genes,

    cd ShotgunUnifrac/CorGE
    pip install .

## Prereqs

  1. Anaconda/miniconda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation)
  2. Snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
  3. (For testing only) PyTest
  4. (For containerization only) Singularity

## Running tests

For the `CorGE` package,

    pytest CorGE/tests

For the tree building Snakemake workflow,

    pytest .tests/
    snakemake --lint
    snakefmt workflow/Snakefile
    
N.B. The linting and formatting steps are optional, you should run them on any contributions you want to make to the codebase, if you want to run snakefmt it must be installed with `pip install snakefmt`

## Running

To download and collect genomes for tree building,

    CorGE collect_genomes ./ --ncbi_species LIST_OF_TXIDS.txt --ncbi_accessions LIST_OF_ACCS.txt --local /path/to/local/db --outgroup 2173

And then to filter out genes of interest and curate everything for tree building,

    CorGE extract_genes /path/to/previous/output --output . --type prot

The default `--type` behavior is 'prot' so that can be left off or switched to 'nucl' if you want to build trees based on nucleotide sequences. Finally to generate the tree, make sure you're in the directory with all the output from the previous step and run,

    snakemake -c --use-conda --conda-prefix .snakemake/

This should output a file called `RAxML_rootedTree.final` which contains the final tree

A worked example is given in the [docs](https://shotgununifrac.readthedocs.io/en/latest/quickstart.html).
