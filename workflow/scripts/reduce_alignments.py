import csv
from io import TextIOWrapper

def read_seqs(f: TextIOWrapper) -> dict:
    d = dict()
    current = ""
    for l in f.readlines():
        if l[0] == ">":
            current = l[2:]
            d[current] = ""
        else:
            d[current] += l.strip()

    return d

#for fp_stats, fp_reduced, fp_seq in zip(snakemake.input.stats, snakemake.output, snakemake.input.seqs):
with open(str(snakemake.input.stats)) as f_stats, open(str(snakemake.output), "w") as f_reduced, open(str(snakemake.input.seqs)) as f_seq:
    r_stats = csv.reader(f_stats, delimiter="\t")
    stats_header = next(r_stats)
    keepers = list()
    for l in r_stats:
        if float(l[2]) < 0.50 and float(l[3]) > 0.0:
            keepers.append((int(l[0]), float(l[3])))

    keepers.sort(key=lambda t: t[1]) # Get 100 lowest entropy columns
    keepers = keepers[:99]
    keepers.sort(key=lambda t: t[0]) # Order by column number to iterate through input fasta
    print(keepers)

    seqs = read_seqs(f_seq)

    for k, v in seqs.items():
        f_reduced.write(f"> {k.strip()}\n")
        for col_num, entropy in keepers:
            f_reduced.write(f"{v[col_num-1]}")
        f_reduced.write("\n")
