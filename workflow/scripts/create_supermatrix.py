quasigenomes = {}
for fp in snakemake.input:
    with open(fp) as f:
        id = ""
        for l in f.readlines():
            if l[0] == ">":
                id = l[2:].strip(" \n\r")
            else:
                try:
                    quasigenomes[id] += l.strip(" \n\r")
                except KeyError:
                    quasigenomes[id] = l.strip(" \n\r")

    max_row_len = max(map(len, quasigenomes.values()))
    for k, v in quasigenomes.items():  # Fill missed genes with blanks
        if len(v) < max_row_len:
            filler = "-" * (max_row_len - len(v))
            quasigenomes[k] += filler

with open(snakemake.output[0], "w") as f:
    for k, v in quasigenomes.items():
        print(f"{len(v)}")
        f.write(f"> {k}\n")
        f.write(f"{v}\n")
