import collections
import csv
import logging
import os
import pyhmmer.plan7
import pyhmmer.easel
import shutil
import sys
import tqdm
import urllib.request
import wget
from io import TextIOWrapper
from warnings import warn

INPUT_FP = ""
OUTPUT_FP = ""
FILTERED_SEQUENCES_FP = ""
MERGED_SEQUENCES_FP = ""

def check_assembly_summary():
    if not os.path.exists(os.path.join(INPUT_FP, "assembly_summary.txt")):
        logging.log(1, "assembly_summary.txt not found, fetching...")
        wget.download("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt", out=INPUT_FP)

        with open(os.path.join(INPUT_FP, "assembly_summary.txt"), "a") as f: # Append 2173 default outgroup
            f.write("GCF_000016525.1\tPRJNA224116\tSAMN02604313\t\trepresentative genome\t420247\t2173\tMethanobrevibacter smithii ATCC 35061\tstrain=ATCC 35061; PS; DSMZ 861\t\tlatest\tComplete Genome\tMajor\tFull\t2007/06/04\tASM1652v1\tWashington University Center for Genome Sciences\tGCA_000016525.1\tidentical\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/525/GCF_000016525.1_ASM1652v1\t\tassembly from type material\tna")

def val_for(acc: str, header: str) -> str:
    check_assembly_summary()
    
    with open(os.path.join(INPUT_FP, "assembly_summary.txt")) as f:
        reader = csv.reader(f, dialect=csv.excel_tab)
        next(reader) # First row is a comment
        headers = next(reader) # This row has the headers
        headers[0] = headers[0][2:]# Remove the "# " from the beginning of the first element

        idIndex = headers.index(header)
        accIndex = headers.index("assembly_accession")

        for line in reader:
            try:
                if line[accIndex] == acc:
                    return line[idIndex]
            except IndexError:
                None # Incomplete assembly_summary entry

    return acc

def txid_for(acc: str) -> str:
    return val_for(acc, "species_taxid")

def strain_for(acc: str) -> str:
    return val_for(acc, "infraspecific_name")

def species_for(acc: str) -> str:
    return val_for(acc, "organism_name")

def write_sequence(out: TextIOWrapper, seqs: TextIOWrapper, query: str, acc: str = None):
    seq = ""
    add = False
    for l in seqs.readlines():
        if add:
            if l[0] != ">":
                seq += l.strip()
            else:
                break
        if query in l and l[0] == '>':
            add = True
            if acc:
                seq += f"> {acc}\n"
            else:
                seq += f"{l.strip()}\n"
    
    out.write(seq + "\n")

# From pyhmmer docs https://pyhmmer.readthedocs.io/en/stable/examples/fetchmgs.html
def run_hmmscan(proteins: list) -> list:
    url = "https://github.com/motu-tool/fetchMGs/raw/master/lib/MG_BitScoreCutoffs.allhits.txt"

    cutoffs = {}
    with urllib.request.urlopen(url) as f:
        for line in csv.reader(TextIOWrapper(f), dialect="excel-tab"):
            if not line[0].startswith("#"):
                cutoffs[line[0]] = float(line[1])

    baseurl = "https://github.com/motu-tool/fetchMGs/raw/master/lib/{}.hmm"

    hmms = []
    for cog in cutoffs:
        with urllib.request.urlopen(baseurl.format(cog)) as f:
            hmm = next(pyhmmer.plan7.HMMFile(f))
            cutoff = cutoffs[hmm.name.decode()]
            hmm.cutoffs.trusted = (cutoff, cutoff)
            hmms.append(hmm)

    Result = collections.namedtuple("Result", ["query", "cog", "bitscore"])
    
    results = []
    for top_hits in pyhmmer.hmmsearch(hmms, proteins, bit_cutoffs="trusted"):
        for hit in top_hits:
            cog = hit.best_domain.alignment.hmm_name.decode()
            results.append(Result(hit.name.decode(), cog, hit.score))

    best_results = {}
    keep_query = set()
    for result in results:
        if result.query in best_results:
            previous_bitscore = best_results[result.query].bitscore
            if result.bitscore > previous_bitscore:
                best_results[result.query] = result
                keep_query.add(result.query)
            elif result.bitscore == previous_bitscore:
                if best_results[result.query].cog != hit.cog:
                    keep_query.remove(result.query)
        else:
            best_results[result.query] = result
            keep_query.add(result.query)

    return [best_results[k] for k in sorted(best_results) if k in keep_query]

def filter_sequences():
    prot_input_fp = os.path.join(INPUT_FP, "protein")
    outgroup_input_fp = os.path.join(INPUT_FP, "outgroup")
    
    prot_list = os.listdir(prot_input_fp)
    prot_fp_list = [os.path.join(prot_input_fp, prot) for prot in prot_list]
    outgroup = os.listdir(outgroup_input_fp)
    for fp in outgroup:
        if ".faa" in fp:
            prot_fp_list.append(os.path.join(outgroup_input_fp, fp))

    logging.log(1, "Filtering sequences...")
    with tqdm.tqdm(total=len(prot_list)) as pbar:
        for prot_fp in prot_fp_list:
            pbar.update(1)
            with pyhmmer.easel.SequenceFile(prot_fp, digital=True) as seqs_file:
                proteins = list(seqs_file)
            
            results = run_hmmscan(proteins)

            acc = prot_fp.split('/')[-1].split('.faa')[0]

            for result in results:
                with open(os.path.join(FILTERED_SEQUENCES_FP, f"{result.cog}__{acc}.faa"), "w") as f:
                    with open(prot_fp) as g:
                        write_sequence(f, g, result.query)

            for result in results[:10]:
                print(result.query, "{:.1f}".format(result.bitscore), result.cog, sep="\t")

def filter_nucl_sequences():
    global FILTERED_SEQUENCES_FP, OUTPUT_FP
    filtered_prot_sequences_fp = FILTERED_SEQUENCES_FP
    FILTERED_SEQUENCES_FP = os.path.join(OUTPUT_FP, "filtered-nucl-sequences")
    try:
        os.mkdir(FILTERED_SEQUENCES_FP)
    except FileExistsError:
        warn("filtered-nucl-sequences output directory already exist and will be overwritten")
        shutil.rmtree(FILTERED_SEQUENCES_FP)
        os.mkdir(FILTERED_SEQUENCES_FP)
    except OSError:
        sys.exit("Problem creating nucl output directories")

    nucl_input_fp = os.path.join(INPUT_FP, "nucleotide")
    outgroup_input_fp = os.path.join(INPUT_FP, "outgroup")
    outgroup_acc = os.listdir(outgroup_input_fp)[0].split(".f")[0]

    for fp in os.listdir(filtered_prot_sequences_fp):
        cog = fp.split("__")[0]
        acc = fp.split("__")[1].split(".faa")[0]
        
        query = ""
        with open(os.path.join(filtered_prot_sequences_fp, fp)) as f:
            query = f.readline().strip()[1:].split(' ')[0]
        
        with open(os.path.join(FILTERED_SEQUENCES_FP, f"{cog}__{acc}.fna"), "w") as f:
            if acc != outgroup_acc:
                with open(os.path.join(nucl_input_fp, f"{acc}.fna")) as g:
                    write_sequence(f, g, query)
            else:
                with open(os.path.join(outgroup_input_fp, f"{acc}.fna")) as g:
                    write_sequence(f, g, query)

def merge_sequences(nt: str):
    all_filtered_seqs = os.listdir(FILTERED_SEQUENCES_FP)
    
    for seq in all_filtered_seqs:
        cog = seq.split("__")[0]
        acc = seq.split("__")[1].split(".f")[0]
        with open(os.path.join(FILTERED_SEQUENCES_FP, seq)) as f:
            with open(os.path.join(MERGED_SEQUENCES_FP, f"{cog}.fasta"), "a") as g:
                if nt == "acc":
                    g.write(f"> {acc}\n")
                elif nt == "txid":
                    g.write(f"> {txid_for(acc)}\n")
                elif nt == "strain":
                    g.write(f"> {strain_for(acc)}\n")
                elif nt == "species":
                    g.write(f"> {species_for(acc)}\n")
                else:
                    sys.exit("Invalid name type supplied")
                g.write(f.readlines()[1])


def write_config(t: str):
    outgroup_input_fp = os.path.join(INPUT_FP, "outgroup")

    all_merged_seqs = os.listdir(MERGED_SEQUENCES_FP)

    cogs = [fp.split(".fasta")[0] for fp in all_merged_seqs]
    
    cfg = "# Config file for tree building pipeline\n# Genes to build trees from\nGENES: ["
    for cog in cogs:
        cfg += f"\"{cog.strip()}\", "
    cfg += "]\n"

    cfg += "# Outgroup to use for rooting, false if outgroup rooting shouldn't be used\n"
    if len(os.listdir(outgroup_input_fp)) == 2:
        
        cfg += f"OUTGROUP: {txid_for(os.listdir(outgroup_input_fp)[0].split('.f')[0])}\n" #TODO: change this to specify outgroup name
    else:
        cfg += "OUTGROUP: false\n"
    
    cfg += f"# File type contained in merged-sequences (prot or nucl)\nTYPE: {t}\n"
    

    with open(os.path.join(OUTPUT_FP, "config.yml"), "w") as f:
        f.write(cfg)

def extract_genes(genomes: str, output: str, output_type: str, name_type: str):
    if not os.path.isdir(genomes):
        sys.exit("Invalid path to collected genomes")
    if not os.path.isdir(output):
        sys.exit("Invalid output path")
    if [fn.split('.fna')[0] for fn in os.listdir(os.path.join(genomes, "nucleotide"))].sort() != \
        [fn.split('.faa')[0] for fn in os.listdir(os.path.join(genomes, "protein"))].sort():
            sys.exit(f"Contents of {os.path.join(genomes, 'nucleotide')} and {os.path.join(genomes, 'protein')} are not perfectly paired")


    global INPUT_FP, OUTPUT_FP, FILTERED_SEQUENCES_FP, MERGED_SEQUENCES_FP
    INPUT_FP = genomes
    OUTPUT_FP = output
    FILTERED_SEQUENCES_FP = os.path.join(output, "filtered-sequences")
    MERGED_SEQUENCES_FP = os.path.join(output, "merged-sequences")

    try:
        os.mkdir(FILTERED_SEQUENCES_FP)
    except FileExistsError:
        warn("filtered-sequences output directory already exist and will be overwritten")
        shutil.rmtree(FILTERED_SEQUENCES_FP)
        os.mkdir(FILTERED_SEQUENCES_FP)
    except OSError:
        sys.exit("Problem creating output directories")

    try:
        os.mkdir(MERGED_SEQUENCES_FP)
    except FileExistsError:
        warn("merged-sequences output directory already exist and will be overwritten")
        shutil.rmtree(MERGED_SEQUENCES_FP)
        os.mkdir(MERGED_SEQUENCES_FP)
    except OSError:
        sys.exit("Problem creating output directories")

    filter_sequences()
    if str(output_type) == 'nucl':
        filter_nucl_sequences()
    merge_sequences(name_type)
    write_config(output_type)
