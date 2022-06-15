import os
import sys
import warnings

def extract_genes(genomes: str, output: str):
    if not os.path.isdir(genomes):
        sys.exit("Invalid path to collected genomes")
    if not os.path.isdir(output):
        sys.exit("Invalid output path")
    
    os.listdir(genomes) # Test read permissions

    filtered_sequences_fp = os.path.join(output, "filtered-sequences")
    merged_sequences_fp = os.path.join(output, "merged-sequences")
    try:
        os.mkdir(filtered_sequences_fp)
        os.mkdir(merged_sequences_fp)
    except FileExistsError:
        warnings.warn("One or more output directories already exist and will be overwritten")
    except OSError:
        sys.exit("Problem creating output directories")