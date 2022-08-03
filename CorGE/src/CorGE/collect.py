import csv
import gzip
import logging
import os
import shutil
import sys
import tqdm
import wget
from io import TextIOWrapper
from warnings import warn

OUTPUT_FP = ""
NUCL_OUTPUT_FP = ""
PROT_OUTPUT_FP = ""
OUTGROUP_OUTPUT_FP = ""

def check_assembly_summary():
    if not os.path.exists(os.path.join(OUTPUT_FP, "assembly_summary.txt")):
        logging.log(1, "assembly_summary.txt not found, fetching...")
        wget.download("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt", out=OUTPUT_FP)

        with open(os.path.join(OUTPUT_FP, "assembly_summary.txt"), "a") as f: # Append 2173 default outgroup
            f.write("GCF_000016525.1\tPRJNA224116\tSAMN02604313\t\trepresentative genome\t420247\t2173\tMethanobrevibacter smithii ATCC 35061\tstrain=ATCC 35061; PS; DSMZ 861\t\tlatest\tComplete Genome\tMajor\tFull\t2007/06/04\tASM1652v1\tWashington University Center for Genome Sciences\tGCA_000016525.1\tidentical\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/525/GCF_000016525.1_ASM1652v1\t\tassembly from type material\tna")


def accession_for(txid: str) -> str or None:
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
                if int(line[idIndex]) == int(txid) and line[refSeqIndex] != "na":
                    return line[accIndex]
            except IndexError:
                None # Incomplete assembly_summary entry

    return None

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
                if line[accIndex] == acc.strip():
                    return line[ftpIndex]
            except IndexError:
                None # Incomplete assembly_summary entry

def download(url_prefix: str, full_acc: str, acc: str, out: str, suffix: str, ext: str):
    try:
        wget.download(f"{url_prefix}/{full_acc}{suffix}{ext}.gz", out=out)
    except Exception as e:
        warn(str(e))
    
    with gzip.open(os.path.join(out, f"{full_acc}{suffix}{ext}.gz"), "rb") as f_in:
        with open(os.path.join(out, f"{acc}{ext}"), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    try:
        os.remove(os.path.join(out, f"{full_acc}{suffix}{ext}.gz"))
    except OSError as e:
        warn(str(e))

def fetch_genome(acc: str or None):
    if acc:
        url_prefix = url_for(acc)
        full_acc = url_prefix.split('/')[-1]
        
        if f"{acc}.fna" in os.listdir(NUCL_OUTPUT_FP):
            logging.log(1, f"Found existing nucleotide file for {acc}")
        else:
            download(url_prefix, full_acc, acc, NUCL_OUTPUT_FP, "_cds_from_genomic", ".fna")

        if f"{acc}.faa" in os.listdir(PROT_OUTPUT_FP):
            logging.log(1, f"Found existing protein file for {acc}")
        else:
            download(url_prefix, full_acc, acc, PROT_OUTPUT_FP, "_protein", ".faa")
    else:
        logging.log(1, "Failed to find accession")

def collect_all(dryrun: bool):
    check_assembly_summary()
    
    accs = []
    ids = []
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
                if int(line[idIndex]) not in ids and line[refSeqIndex] != "na":
                    ids.append(int(line[idIndex]))
                    accs.append(line[accIndex])
            except IndexError:
                None # Incomplete assembly_summary entry
    
    if dryrun:
        print(f"Accessions for all RefSeq species to be collected: {len(accs)}")
        #print(", ".join(accs))
    else:
        for acc in accs:
            fetch_genome(acc)

def collect_ncbi_species(ids: TextIOWrapper, dryrun: bool):
    accs = []
    for id in ids.readlines():
        if dryrun:
            accs.append(accession_for(id))
        else:
            fetch_genome(accession_for(id))

    if dryrun:
        accs = list(filter(None, accs))
        print(f"Accessions based on input taxon ids to be collected: {len(accs)}")
        print(", ".join(accs))

def collect_ncbi_accessions(accs: TextIOWrapper, dryrun: bool):
    accs = []
    for acc in accs.readlines():
        if dryrun:
            accs.append(acc)
        else:
            fetch_genome(acc)

    if dryrun:
        print(f"Accessions based on input accessions to be collected: {len(accs)}")
        print(", ".join(accs))

def collect_local(dir: str, dryrun: bool):
    nucl_list = list()
    prot_list = list()
    for fp in os.listdir(dir):
        if ".fna" in fp:
            nucl_list.append(fp)
        if ".faa" in fp:
            prot_list.append(fp)

    fps = list()
    for n in nucl_list:
        if n.replace(".fna", ".faa") in prot_list:
            if dryrun:
                fps.append(n.replace(".fna", ""))
            else:
                shutil.copyfile(os.path.join(dir, n), os.path.join(NUCL_OUTPUT_FP, n))
                shutil.copyfile(os.path.join(dir, n.replace(".fna", ".faa")), os.path.join(PROT_OUTPUT_FP, n.replace(".fna", ".faa")))

    if dryrun:
        print(f"Local genomes to be collected: {len(fps)}")
        print(", ".join(fps))

def collect_outgroup(id: str, dryrun: bool):
    if dryrun:
        print("Outgroup genome to be collected:")
        print(accession_for(id))
    else:
        if len(os.listdir(OUTGROUP_OUTPUT_FP)):
            warn("Overwriting existing outgroup files")
            shutil.rmtree(OUTGROUP_OUTPUT_FP)
            os.mkdir(OUTGROUP_OUTPUT_FP)
        if id.isdigit():
            id = accession_for(id)
            fetch_genome(id)
        
        try:
            os.rename(os.path.join(NUCL_OUTPUT_FP, f"{id}.fna"), os.path.join(OUTGROUP_OUTPUT_FP, f"{id}.fna"))
            os.rename(os.path.join(PROT_OUTPUT_FP, f"{id}.faa"), os.path.join(OUTGROUP_OUTPUT_FP, f"{id}.faa"))
        except OSError:
            sys.exit(f"Couldn't find genome files corresponding to outgroup {id}")

def collect_genomes(output: str, all: bool, ncbi_species: TextIOWrapper or None, ncbi_accessions: TextIOWrapper or None, local: str, outgroup: str, dryrun: bool):
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
    None if os.path.isdir(NUCL_OUTPUT_FP) else os.mkdir(NUCL_OUTPUT_FP)
    None if os.path.isdir(PROT_OUTPUT_FP) else os.mkdir(PROT_OUTPUT_FP)
    None if os.path.isdir(OUTGROUP_OUTPUT_FP) else os.mkdir(OUTGROUP_OUTPUT_FP)

    if all:
        collect_all(dryrun)
    if ncbi_species:
        collect_ncbi_species(ncbi_species, dryrun)
    if ncbi_accessions:
        collect_ncbi_accessions(ncbi_accessions, dryrun)
    if local != "":
        collect_local(local, dryrun)
    if outgroup != "None":
        collect_outgroup(outgroup, dryrun)