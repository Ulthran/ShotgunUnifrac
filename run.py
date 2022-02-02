#!~/anaconda3/bin/ python3
# Wrapper for the tree building pipeline
# @param input is the Kraken input file

import csv
import sys
import os

# Write config file for this run
input = sys.argv[1]
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

os.system("snakemake -c --use-conda")
