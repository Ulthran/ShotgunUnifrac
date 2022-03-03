# Downloads genome from NCBI for a given genome id
import csv
import gzip
from io import TextIOWrapper
import os
import shutil
import tqdm
import wget
from typing import NoReturn, Tuple

# Checks for the existence of data/assembly_summary.txt
# @return is True if the file exists, False otherwise
def check_for_assembly() -> bool:
    return os.path.exists("data/assembly_summary.txt")

# Downloads assembly_summary.txt if it doens't already exist
# @return is the file downloaded by wget
def download_assembly() -> int:
    print("data/assembly_summary.txt not found, fetching...")
    return wget.download("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt", out="data/")

# Creates a file of the best genome accessions for a given list of species-level taxon ids
# @param txids is the list of species-level taxon ids
# @param naTol is the 'na'-tolerance, if True the function will include accessions with refseq_category of na, False by default
# @return is the list of txids, where each one that a genome accession was found for has been removed
def find_genome_accessions(txids: list, naTol: bool = False) -> list:
    with tqdm.tqdm(total=os.path.getsize("data/assembly_summary.txt")) as pbar: # Add progress bar
        with open("data/run_assembly.txt", "a") as run_assembly:
            with open("data/assembly_summary.txt") as assembly:
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

# Prepares the run_assembly.txt file for downloading genomes
# @param inputFile is the path to the list of taxon ids
# @param outputDir is the output location for run_assembly.txt
# @param logF is the log file to write to
# @return is the list of taxon ids for which a reference/representative genome couldn't be found
def prepare_run_assembly(inputFile: str, outputDir: str, logF: TextIOWrapper) -> list:
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
    with open(os.path.join(outputDir, "data/run_assembly.txt"), "w") as run_assembly:
        with open("data/assembly_summary.txt") as assembly:
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
   
    return txids

# Downloads genomes from NCBI for all the accessions in run_assembly.txt
# @param outputDir is the location to download genomes
# @param logF is the log file to write to
# @return is a tuple containing the ids of successfully downloaded genomes and those of failed genomes
def download_genomes(outputDir: str, logF: TextIOWrapper) -> Tuple[list, list]:
    downloaded_genome_ids = []
    failed_genome_ids = []
    ncbi_dir = os.path.join(outputDir, "output/ncbi/")
    run_assembly_path = os.path.join(outputDir, "data/run_assembly.txt")
    with open (run_assembly_path) as run_assembly:
        run_assembly_reader = csv.reader(run_assembly, dialect=csv.excel_tab)
        firstLine = next(run_assembly_reader)
        ftpIndex = firstLine.index("ftp_path")

        for line in run_assembly_reader:
            genomeId = line[ftpIndex].split("/")[-1]
            for ext in ["_cds_from_genomic.fna.gz", "_rna_from_genomic.fna.gz"]:
                if not os.path.isfile(ncbi_dir + genomeId + ext[:-6] + "fasta"):
                    url = line[ftpIndex] + "/" + genomeId + ext
                    wget.download(url, out=ncbi_dir)
                    with gzip.open(ncbi_dir + genomeId + ext, "rb") as f_in:
                        with open(ncbi_dir + genomeId + ext[:-6] + "fasta", "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    try:
                        os.remove(ncbi_dir + genomeId + ext)
                    except OSError as e:
                        print(e)
                else:
                    logF.write("Found " + ncbi_dir + genomeId + ext[:-6] + "fasta, skipping download\n")
            if os.path.isfile(ncbi_dir + genomeId + "_cds_from_genomic.fasta") and os.path.isfile(ncbi_dir + genomeId + "_rna_from_genomic.fasta"):
                downloaded_genome_ids.append(genomeId)
            else:
                failed_genome_ids.append(genomeId)
    return downloaded_genome_ids, failed_genome_ids

# Wrapper for download_genes from snakemake
# @param output is the output file containing the genome id
# @return is undefined
def download_genes_sm(output: str) -> NoReturn:
    outputs = output.split(" ")
    for out in outputs:
        genomeId = out[5:-23]
        ext = out[-23:-5] + "fna.gz"
        download_genome(genomeId, ext)













