###
# All code here taken from okfasta
# https://github.com/kylebittinger/okfasta
###
from io import StringIO
from .msa import MSA

def parse_fasta(f, trim_desc=False):
    f = iter(f)
    desc = next(f).strip()[1:]
    if trim_desc:
        desc = desc.split()[0]
    seq = StringIO()
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            yield desc, seq.getvalue()
            desc = line[1:]
            if trim_desc:
                desc = desc.split()[0]
            seq = StringIO()
        else:
            line = line.replace(" ", "").replace("U", "T").replace(".", "-")
            seq.write(line)
    yield desc, seq.getvalue()
with open(str(snakemake.input)) as f_in, open(str(snakemake.output), "w") as f_out:
    seqs = parse_fasta(f_in)
    msa = MSA.from_seqs(seqs)
    f_out.write("\t".join(msa.column_stats_header))
    f_out.write("\n")
    for stats_result in msa.column_stats():
        stats_values = stats_result.values()
        f_out.write(msa.column_stats_fmt.format(*stats_values))
        f_out.write("\n")