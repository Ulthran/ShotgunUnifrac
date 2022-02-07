#!~/anaconda3/bin/ python3
# Wrapper for the tree building pipeline
# @param input is the Kraken input file
# @param test is a boolean to determine whether or not to run tests before
# running the pipeline, defaults to False
# @param singularity is a boolean to determine whether or not to use the
# singularity image when running the pipeline, defaults to False

import csv
import sys
import os
import shutil
import pytest

input = sys.argv[1]
test = False
singularity = False
try:
    test = sys.argv[2]
except IndexError:
    None
try:
    singularity = sys.argv[3]
except IndexError:
    None

# Write config file for this run
with open(input, newline="") as inputF:
    tsv = csv.reader(inputF, dialect=csv.excel_tab)
    idIndex = next(tsv).index("ftp_path")
    with open("config.yml", "w") as config:
        config.write("# Config file for tree building pipeline\n")
        config.write("# Kraken output file\n")
        config.write("KRAKEN: \"" + input + "\"\n")
        config.write("# Genome accessions (NCBI)\n")
        config.write("IDS: [")
        for line in tsv:
            val = line[idIndex].split("/")[-1]
            config.write("\"" + val + "\", ")
        config.write("]")
print("Created config.yml")

# Run tests
if(test):
    retcode = pytest.main(["-x", ".tests/"])
    os.system("snakemake --lint")

# Check for singularity and install if not there, then run with singularity
if(singularity):
    if(not shutil.which(singularity)):
        os.system("yes | pip install singularity")
    os.system("snakemake -c --use-conda --use-singularity")
else:
    os.system("snakemake -c --use-conda")
