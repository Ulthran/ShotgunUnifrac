# ShotgunUnifrac

## Install

    git clone git@github.com:Ulthran/ShotgunUnifrac.git

## Prereqs

  1. Anaconda/miniconda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation)
  2. Snakemake
  3. (For testing only) PyTest

## Running tests

    pytest .tests/
    snakemake --lint

## Running

To generate tree:

    python run.py workflow/data/kraken2_standard_20200204_select_genomes.txt True True

Arguments are:
  1. str: Path to kraken report
  2. bool[False]: Run tests before pipeline
  3. bool[False]: Run pipeline with singularity