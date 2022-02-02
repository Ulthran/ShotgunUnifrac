# ShotgunUnifrac

## Install

    git clone git@github.com:Ulthran/ShotgunUnifrac.git
    cd ShotgunUnifrac
    conda env create -f environment.yml

## Running

To generate tree:

    conda activate shotguntree_env
    python run.py data/kraken2_standard_20200204_select_genomes.txt

## Running tests

    pytest .tests/