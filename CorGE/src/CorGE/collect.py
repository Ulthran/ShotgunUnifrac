import csv
import gzip
import os
import shutil
import sys
import wget
from io import TextIOWrapper
from logging import log
from warnings import warn

OUTPUT_FP = ""
NUCL_OUTPUT_FP = ""
PROT_OUTPUT_FP = ""
OUTGROUP_OUTPUT_FP = ""

def check_assembly_summary():
    if not os.path.exists(os.path.join(OUTPUT_FP, "assembly_summary.txt")):
        log("assembly_summary.txt not found, fetching...")
        wget.download("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt", out=OUTPUT_FP)

def accession_for(txid: str) -> str:
    check_assembly_summary()
    
    with open(os.path.join(OUTPUT_FP, "assembly_summary.txt")) as f:
        reader = csv.reader(f, dialect=csv.excel_tab)
        next(reader) # First row is a comment
        headers = next(reader) # This row has the headers
        headers[0] = headers[0][2:]# Remove the "# " from the beginning of the first element

        idIndex = headers.index("species_taxid")
        accIndex = headers.index("assembly_accession")
        refSeqIndex = headers.index("refseq_category")

        for line in reader:
            try:
                if line[idIndex] == txid and (line[refSeqIndex] == "reference genome" or line[refSeqIndex] == "representative genome"):
                    return line[accIndex]
            except IndexError:
                None # Incomplete assembly_summary entry

def url_for(acc: str) -> str:
    check_assembly_summary()

    with open(os.path.join(OUTPUT_FP, "assembly_summary.txt")) as f:
        reader = csv.reader(f, dialect=csv.excel_tab)
        next(reader) # First row is a comment
        headers = next(reader) # This row has the headers
        headers[0] = headers[0][2:]# Remove the "# " from the beginning of the first element

        accIndex = headers.index("assembly_accession")
        ftpIndex = headers.index("ftp_path")

        for line in reader:
            try:
                if line[accIndex] == acc:
                    return line[ftpIndex]
            except IndexError:
                None # Incomplete assembly_summary entry

def download(url_prefix: str, full_acc: str, acc: str, out: str, suffix: str, ext: str):
    try:
        wget.download(f"{url_prefix}/{full_acc}{suffix}{ext}.gz", out=out)
    except Exception as e:
        warn(str(e))
    
    with gzip.open(os.path.join(out, f"{full_acc}{suffix}{ext}.gz"), "rb") as f_in:
        with open(os.path.join(out, f"{acc}.fna"), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    try:
        os.remove(os.path.join(out, f"{full_acc}{suffix}{ext}.gz"))
    except OSError as e:
        warn(str(e))

def fetch_genome(acc: str):
    url_prefix = url_for(acc)
    full_acc = url_prefix.split('/')[-1]
    
    if f"{acc}.fna" in os.listdir(NUCL_OUTPUT_FP):
        log(f"Found existing nucleotide file for {acc}")
    else:
        download(url_prefix, full_acc, acc, NUCL_OUTPUT_FP, "_cds_from_genomic", ".fna")

    if f"{acc}.faa" in os.listdir(PROT_OUTPUT_FP):
        log(f"Found existing protein file for {acc}")
    else:
        download(url_prefix, full_acc, acc, PROT_OUTPUT_FP, "_protein", ".faa")

def collect_ncbi_species(ids: TextIOWrapper):
    for id in ids.readlines():
        fetch_genome(accession_for(id))

def collect_ncbi_accessions(accs: TextIOWrapper):
    for acc in accs.readlines():
        fetch_genome(acc)

def collect_local(dir: str):
    nucl_list = list()
    prot_list = list()
    for fp in os.listdir(dir):
        if ".fna" in fp:
            nucl_list.append(fp)
        if ".faa" in fp:
            prot_list.append(fp)

    for n in nucl_list:
        if n in prot_list:
            shutil.copyfile(os.path.join(dir, f"{n}.fna"), os.path.join(NUCL_OUTPUT_FP, f"{n}.fna"))
            shutil.copyfile(os.path.join(dir, f"{n}.faa"), os.path.join(PROT_OUTPUT_FP, f"{n}.faa"))

def collect_outgroup(id: str):
    if len(os.listdir(OUTGROUP_OUTPUT_FP)):
        warn("Overwriting existing outgroup files")
        shutil.rmtree(OUTGROUP_OUTPUT_FP)
        os.mkdir(OUTGROUP_OUTPUT_FP)
    if id.isdigit():
        fetch_genome(accession_for(id))
    
    try:
        os.rename(os.path.join(NUCL_OUTPUT_FP, f"{id}.fna"), os.path.join(OUTGROUP_OUTPUT_FP, f"{id}.fna"))
        os.rename(os.path.join(PROT_OUTPUT_FP, f"{id}.faa"), os.path.join(OUTGROUP_OUTPUT_FP, f"{id}.faa"))
    except OSError:
        sys.exit(f"Couldn't find genome files corresponding to outgroup {id}")

def collect_genomes(output: str, ncbi_species: TextIOWrapper or None, ncbi_accessions: TextIOWrapper or None, local: str, outgroup: str):
    if not os.path.isdir(output):
        try:
            os.makedirs(output)
        except OSError:
            sys.exit("Could not create outputs in specified directory")
    
    global OUTPUT_FP, NUCL_OUTPUT_FP, PROT_OUTPUT_FP, OUTGROUP_OUTPUT_FP
    OUTPUT_FP = output
    NUCL_OUTPUT_FP = os.path.join(output, "nucleotide")
    PROT_OUTPUT_FP = os.path.join(output, "protein")
    OUTGROUP_OUTPUT_FP = os.path.join(output, "outgroup")
    None if not os.path.isdir(NUCL_OUTPUT_FP) else os.mkdir(NUCL_OUTPUT_FP)
    None if not os.path.isdir(PROT_OUTPUT_FP) else os.mkdir(PROT_OUTPUT_FP)
    None if not os.path.isdir(OUTGROUP_OUTPUT_FP) else os.mkdir(OUTGROUP_OUTPUT_FP)

    if ncbi_species:
        collect_ncbi_species(ncbi_species)
    if ncbi_accessions:
        collect_ncbi_accessions(ncbi_accessions)
    if local != "":
        collect_local(local)
    if outgroup != "None":
        collect_outgroup(outgroup)