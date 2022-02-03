# ShotgunUnifrac

## Install

    git clone git@github.com:Ulthran/ShotgunUnifrac.git

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