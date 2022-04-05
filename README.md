# ShotgunUnifrac

A dual use program for downloading and extracting genes from NCBI and for creating phylogenetic trees for many marker genes and merging the results into one

## Install

    git clone git@github.com:Ulthran/ShotgunUnifrac.git

To install the `genetools` library for downloading, extracting, and merging genes from NCBI,

    cd ShotgunUnifrac/workflow/scripts
    pip install .

## Prereqs

  1. Anaconda/miniconda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation)
  2. Snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
  3. (For testing only) PyTest
  4. (For containerization only) Singularity

## Running tests

For the `genetools` package,

    pytest workflow/scripts/tests

For the tree building Snakemake workflow,

    pytest .tests/
    snakemake --lint
    snakefmt workflow/Snakefile
    
N.B. The linting and formatting steps are optional, you should run them on any contributions you want to make to the codebase, if you want to run snakefmt it must be installed with `pip install snakefmt`

## Running

To download and curate data for tree building,

    cd workflow/scripts
    genetools /path/to/taxonid/list -f /path/to/gene/list -m --write_config

N.B. You can include the `--remove_temp` flag as well to remove all the files except the final output, in this case the merged gene files

This should download reference/representative genomes from NCBI for each taxon, extract each gene from each of those genomes into its own file (-f), and then combine those into one file per gene (-m) and write a config for the tree building pipeline (--write_config)

To generate tree:

    cd ../..
    cp workflow/scripts/output/config.yml .
    snakemake --use-conda

This should output a file called `RAxML_rootedTree.final` which contains the final tree

A worked example is given in the [wiki](https://github.com/Ulthran/ShotgunUnifrac/wiki)
