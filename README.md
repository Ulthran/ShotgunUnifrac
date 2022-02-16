# ShotgunUnifrac

## Install

    git clone git@github.com:Ulthran/ShotgunUnifrac.git

## Prereqs

  1. Anaconda/miniconda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation)
  2. Snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
  3. (For testing only) PyTest
  4. (For containerization only) Singularity

N.B. PyTest and Singularity should be installed by `run.py` if they aren't detected on your system

## Running tests

    pytest .tests/
    snakemake --lint
    snakefmt workflow/Snakefile
    
N.B. The linting and formatting steps are optional, you should run them on any contributions you want to make to the codebase, if you want to run snakefmt it must be installed with `pip install snakefmt`

## Running

To generate tree:

    python run.py workflow/data/TEST_TXIDS.txt workflow/data/TEST_GENES.txt True True

Arguments are:
  1. str: Path to list of taxon ids (one per line)
  2. str: Path to list of genes (one per line)
  3. bool[False]: Run tests before pipeline
  4. bool[False]: Run pipeline with singularity (containerized)

A worked example is given in the [wiki](https://github.com/Ulthran/ShotgunUnifrac/wiki)
