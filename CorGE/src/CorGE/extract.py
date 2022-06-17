import collections
import csv
import imp
import logging
import os
from CorGE.src.CorGE.collect import NUCL_OUTPUT_FP
import pyhmmer.plan7
import pyhmmer.easel
import shutil
import sys
import tqdm
import urllib.request
from io import TextIOWrapper
from warnings import warn
#from .command import FileType

from enum import Enum
class FileType(Enum):
    prot = 'prot'
    nucl = 'nucl'

    def __str__(self):
        return self.value

OUTPUT_FP = ""
FILTERED_SEQUENCES_FP = ""
MERGED_SEQUENCES_FP = ""

def write_sequence(out: TextIOWrapper, seqs: TextIOWrapper, query: str, acc: str = None):
    seq = ""
    add = False
    for l in seqs.readlines():
        if add:
            if l[0] != ">":
                seq += l.strip()
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

def filter_sequences(dir: str):
    prot_input_fp = os.path.join(dir, "protein")
    outgroup_input_fp = os.path.join(dir, "outgroup")
    
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

def filter_nucl_sequences(dir: str):
    global FILTERED_SEQUENCES_FP
    filtered_prot_sequences_fp = FILTERED_SEQUENCES_FP
    FILTERED_SEQUENCES_FP = os.path.join(OUTPUT_FP, "filtered-nucl-sequences")
    nucl_input_fp = os.path.join(dir, "nucleotide")

    for fp in os.listdir(filtered_prot_sequences_fp):
        cog = fp.split("__")[0]
        acc = fp.split("__")[1].split(".faa")[0]
        
        query = ""
        with open(os.path.join(filtered_prot_sequences_fp, fp)) as f:
            query = f.readline().strip()[1:].split(' ')[0]
        
        with open(os.path.join(FILTERED_SEQUENCES_FP, f"{cog}__{acc}.fna"), "w") as f:
            with open(os.path.join(nucl_input_fp, f"{acc}.fna")) as g:
                write_sequence(f, g, query)

def merge_sequences():
    all_filtered_seqs = os.listdir(FILTERED_SEQUENCES_FP)
    
    for seq in all_filtered_seqs:
        cog = seq.split("__")[0]
        acc = seq.split("__")[1]
        with open(os.path.join(FILTERED_SEQUENCES_FP, seq)) as f:
            with open(os.path.join(MERGED_SEQUENCES_FP, f"{cog}.fasta"), "a") as g:
                g.write(f"> {acc}\n")
                g.write(f.readlines()[1])


def write_config(dir: str, t: FileType):
    OUTGROUP_INPUT_FP = os.path.join(dir, "outgroup")

    all_merged_seqs = os.listdir(MERGED_SEQUENCES_FP)

    cogs = [fp.split(".faa")[0] for fp in all_merged_seqs]
    
    cfg = "# Config file for tree building pipeline\n# Genes to build trees from\nGENES: ["
    for cog in cogs:
        cfg += f"\"{cog.strip()}\", "
    cfg += "]\n"

    cfg += "# Outgroup to use for rooting, false if outgroup rooting shouldn't be used\n"
    if len(os.listdir(OUTGROUP_INPUT_FP)) == 2:
        cfg += "OUTGROUP: true" #TODO: change this to specify outgroup name
    else:
        cfg += "OUTGROUP: false"
    
    cfg += f"# File type contained in merged-sequences (prot or nucl)\nTYPE: {str(t)}"
    

    with open(os.path.join(OUTPUT_FP, "config.yml"), "w") as f:
        f.write(cfg)

def extract_genes(genomes: str, output: str, output_type: FileType):
    if not os.path.isdir(genomes):
        sys.exit("Invalid path to collected genomes")
    if not os.path.isdir(output):
        sys.exit("Invalid output path")
    if [fn.split('.fna')[0] for fn in os.listdir(os.path.join(genomes, "nucleotide"))].sort() != \
        [fn.split('.faa')[0] for fn in os.listdir(os.path.join(genomes, "protein"))].sort():
            sys.exit(f"Contents of {os.path.join(genomes, 'nucleotide')} and {os.path.join(genomes, 'protein')} are not perfectly paired")


    global OUTPUT_FP, FILTERED_SEQUENCES_FP, MERGED_SEQUENCES_FP
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

    filter_sequences(genomes)
    if str(output_type) == 'nucl':
        filter_nucl_sequences(genomes)
    merge_sequences()
    write_config(genomes, output_type)
