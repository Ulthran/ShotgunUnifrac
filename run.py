#!~/anaconda3/bin/ python3
# Wrapper for the tree building pipeline
# @param inputFile is the list of species-level taxon ids (one id per line)
# @param geneFile is the list of genes to search for (one per line)
# @param test is a boolean to determine whether or not to run tests before
# running the pipeline, defaults to False
# @param singularity is a boolean to determine whether or not to use the
# singularity image when running the pipeline, defaults to False

import csv
import sys
import os
import shutil
from typing import Tuple
try:
    import tqdm
except ImportError:
    os.system("pip install tqdm")
    import tqdm
import workflow.scripts.downloadGenes as dg
from collections import Counter

### Parse command line arguments

inputFile: str = sys.argv[1]
geneFile: str = sys.argv[2]
test: bool = False
singularity: bool = False
try:
    test = sys.argv[3]
except IndexError:
    None
try:
    singularity = sys.argv[4]
except IndexError:
    None

logF = open("logs/main.log", "w")
print("Writing logs to: " + logF.name)

print("Taxon ID list file path: " + inputFile)
print("Gene list file path: " + geneFile)
print("Run tests: " + str(test))
print("Run containerized: " + str(singularity))
logF.write("Taxon ID list file path: " + inputFile + "\n")
logF.write("Gene list file path: " + geneFile + "\n")
logF.write("Run tests: " + str(test) + "\n")
logF.write("Run containerized: " + str(singularity) + "\n")

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
                            if line[idIndex] == txid and (line[refSeqIndex] != "na" or (line[lvlIndex] == "Complete Genome" and naTol)):
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

# Create a list of txids from file of tab delimited columns with txid in the first position
# Ignores lines starting with '#'
txids = []
with open(inputFile) as inputF:
    data = csv.reader(inputF, dialect=csv.excel_tab)
    for line in data:
        if line[0][0] != "#":
            txids.append(line[0])

# Overwrite existing file with header line
with open("run_assembly.txt", "w") as run_assembly:
    with open("workflow/data/assembly_summary.txt") as assembly:
        reader = csv.reader(assembly, dialect=csv.excel_tab)
        next(reader) # First row is a random comment
        firstLine = next(reader) # This row has the headers
        firstLine[0] = firstLine[0][2:]# Remove the "# " from the beginning of the first element

        run_assembly.write("\t".join(firstLine) + "\n") # Write the first line of run_assembly.txt with column identifiers

txids = find_genome_accessions(txids)

if txids:
    print("Could not find suitable entries for:\n")
    logF.write("Could not find suitable entries for:\n")
    for txid in txids:
        print(txid)
        logF.write(txid + "\n")

### Download and curate data for run

None if os.path.isdir("ncbi/") else os.system("mkdir ncbi/")
None if os.path.isdir("sequences/") else os.system("mkdir sequences/")
None if os.path.isdir("merged-sequences/") else os.system("mkdir merged-sequences/")

# Download genomes
downloaded_genome_ids, failed_genome_ids = dg.download_genomes()

if failed_genome_ids:
    logF.write("Failed to download genome(s) for:\n")
    for failure in failed_genome_ids:
        logF.write(failure + "\n")

# Extract genes
gene_counter = dg.extract_genes(geneFile, downloaded_genome_ids)

logF.write("Gene extraction successes (total=" + str(len(downloaded_genome_ids)) + "):\n")
for gene in gene_counter.most_common():
    logF.write(gene[0] + ": " + str(gene[1]) + "\n")
logF.write("Ignoring genes with fewer than 4 sequences...")

# Merge genes
gene_list = dg.merge_genes(gene_counter)

# Clean up
os.system("rm -r sequences/")

### Write config file for this run

print("Writing config.yml for this run")
with open("config.yml", "w") as config:
    config.write("# Config file for tree building pipeline\n")
    config.write("# Genes to build trees from\n")
    config.write("GENES: [")
    for gene in gene_list:
        config.write("\"" + gene + "\", ")
    config.write("]")
print("Created config.yml")

### Check for pytest and install if not there, then run tests

if(test):
    if(not shutil.which(pytest)):
        print("PyTest not installed, installing...")
        os.system("yes | pip install -U pytest")
    import pytest
    retcode = pytest.main(["-x", ".tests/"])
    os.system("snakemake --lint")

### Check for singularity and install if not there, then run with singularity, else run without

if(singularity):
    if(not shutil.which(singularity)):
        print("Singularity not installed, installing...")
        os.system("yes | pip install singularity")
    os.system("snakemake -c --use-conda --use-singularity")
else:
    None
    os.system("snakemake -c --use-conda")

### Post-run summary

logF.close()