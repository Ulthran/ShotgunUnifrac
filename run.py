#!~/anaconda3/bin/ python3
# Wrapper for the tree building pipeline
# @param input is the list of species-level taxon ids (one id per line)
# @param test is a boolean to determine whether or not to run tests before
# running the pipeline, defaults to False
# @param singularity is a boolean to determine whether or not to use the
# singularity image when running the pipeline, defaults to False

import csv
import sys
import os
import shutil
import pytest
import tqdm 

### Parse command line arguments
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

# Checks for the existence of workflow/data/assembly_summary.txt
# @return is True if the file exists, False otherwise
def check_for_assembly() -> bool:
    return os.path.exists("workflow/data/assembly_summary.txt")

# Downloads assembly_summary.txt if it doens't already exist
# @return is the exit status of os.system call to wget
def download_assembly() -> int:
    print("workflow/data/assembly_summary.txt not found, fetching...")
    return os.system("wget -P workflow/data/ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt")

# Creates a file of the best genome accessions for a given list of species-level taxon ids
# @param txids is the list of species-level taxon ids
# @param naTol is the 'na'-tolerance, if True the function will include accessions with refseq_category of na, False by default
# @return is the list of txids, where each one that a genome accession was found for has been removed
def find_genome_accessions(txids: list, naTol: bool = False) -> list:
    with tqdm.tqdm(total=os.path.getsize("workflow/data/assembly_summary.txt")) as pbar: # Add progress bar
        with open("run_assembly.txt", "a") as run_assembly:
            with open("workflow/data/assembly_summary.txt") as assembly:
                reader = csv.reader(assembly, dialect=csv.excel_tab)
                next(reader) # First row is a random comment
                firstLine = next(reader) # This row has the headers
                firstLine[0] = firstLine[0][2:]# Remove the "# " from the beginning of the first element

                run_assembly.write("\t".join(firstLine) + "\n") # Write the first line of run_assembly.txt with column identifiers

                idIndex = firstLine.index("species_taxid")
                accIndex = firstLine.index("assembly_accession")
                lvlIndex = firstLine.index("assembly_level")
                refSeqIndex = firstLine.index("refseq_category")

                # Write each line of assembly_summary.txt that has a genome assembly we want to run_assembly.txt
                for line in reader:
                    pbar.update(len(line)*14) # Multiplying by 14 cause that's what makes the progress bar accurate, issue with line coming from
                                            # a csv.reader object instead of a file object
                    for txid in txids:
                        try:
                            if line[idIndex] == txid[0] and (line[refSeqIndex] != "na" or (line[lvlIndex] == "Complete Genome" and naTol)):
                                run_assembly.write("\t".join(line) + "\n")
                                txids.remove(txid)
                                break
                        except IndexError:
                            None # Incomplete entry in assembly_summary.txt
    return txids

# Check for existence of assembly_summary.txt
None if check_for_assembly() else download_assembly()

### Write run_assembly.txt for this run
print("Writing run_assembly.txt for this run, this could take a minute")

# Create a list of txids from file
txids = []
with open(input) as inputF:
    data = csv.reader(inputF)
    for txid in data:
        txids.append(txid)

with open("run_assembly.txt", "w") as run_assembly:
    run_assembly.write("") # Overwrite existing file, find_genome_accessions works in append mode

txids = find_genome_accessions(txids)

if txids:
    print("Could not find suitable entries for:\n")
    for txid in txids:
        print(txid[0])
    print("Trying with 'na'-tolerance...")
    txids = find_genome_accessions(txids, True)

if txids:
    print("Could not find suitable entries for:\n")
    for txid in txids:
        print(txid[0])
    print("\nNothing left to try. Quitting...")
    quit()

### Write config file for this run
print("Writing config.yml for this run")
with open("run_assembly.txt") as run_assembly:
    tsv = csv.reader(run_assembly, dialect=csv.excel_tab)
    
    genomeIdIndex = next(tsv).index("ftp_path")
    with open("config.yml", "w") as config:
        config.write("# Config file for tree building pipeline\n")
        config.write("# Genome accessions (NCBI)\n")
        config.write("IDS: [")
        for line in tsv:
            val = line[genomeIdIndex].split("/")[-1]
            config.write("\"" + val + "\", ")
        config.write("]\n")
        config.write("# Genes to build trees from\n")
        config.write("GENES: [\"secE\", \"secG\", \"secY\", \"smpB\", \"tsaE\", \"yajC\"]") # Hardcoded for now, should read from data file
print("Created config.yml")

### Run tests
if(test):
    retcode = pytest.main(["-x", ".tests/"])
    os.system("snakemake --lint")

### Check for singularity and install if not there, then run with singularity, else run without
if(singularity):
    if(not shutil.which(singularity)):
        os.system("yes | pip install singularity")
    os.system("snakemake -c --use-conda --use-singularity")
else:
    None
    os.system("snakemake -c --use-conda")
