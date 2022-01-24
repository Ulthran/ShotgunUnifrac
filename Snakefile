SPECIES = [
    "Bifidobacterium longum",
    "Bacteroides dorei",
    "Bacteroides vulgatus",
    "Prevotella copri",
    "Clostridium scindens",
    "Faecalibacterium prausnitzii",
    "Eubacterium rectale",
    "Veillonella parvula",
    "Escherichia coli",
    "Klebsiella pneumoniae"
]
IDS = [
    "GCF_000196555.1_ASM19655v1",
    "GCF_902387545.1_UHGG_MGYG-HGUT-02478",
    "GCF_020885855.1_ASM2088585v1",
    "GCF_020735445.1_ASM2073544v1",
    "GCF_020892115.1_ASM2089211v1",
    "GCF_003312465.1_ASM331246v1",
    
    "GCF_900186885.1_48903_D01",
    "GCF_000005845.2_ASM584v2",
    "GCF_000240185.1_ASM24018v2"
]
GENES = ["secE"]

rule all:
    input: "RAxML_bestTree.rooted"

rule all_gene_files:
    input: expand("sequences/{gene}.fasta", gene=GENES)

rule all_filtered_files:
    input: expand("filtered-sequences/{gene}.fasta", gene=GENES)

rule get_ncbi_sequences:
    input:
        "data/sample.fastq"
    output:
        temp("ncbi/{id}_cds_from_genomic.fasta")
    shell:
        "./download-genes.sh {output}"

rule extract_marker_genes:
    input:
        "ncbi/{id}_cds_from_genomic.fasta"
    output:
        temp("sequences/{gene}__{id}.fasta")
    shell:
        "python3 -c 'import filterGenes; filterGenes.filter_seq_genes_sm(\"ncbi/{wildcards.id}_cds_from_genomic.fasta\", \"{output}\")' > {output}"

rule merge_marker_genes:
    input:
        expand("sequences/{gene}__{id}.fasta", id=IDS, gene=GENES)
    output:
        "sequences/{gene}.fasta"
    shell:
        "cat sequences/{wildcards.gene}__*.fasta > {output}"

rule align_fasta:
    input:
        "sequences/{gene}.fasta"
    output:
        temp("aligned-sequences/{gene}.fasta")
    shell:
        "muscle -in sequences/{wildcards.gene}.fasta -out aligned-sequences/{wildcards.gene}.fasta"

rule filter_columns:
    # Identity op for now, okfasta
    input:
        "aligned-sequences/{gene}.fasta"
    output:
        "filtered-sequences/{gene}.fasta"
    shell:
        "cp aligned-sequences/{wildcards.gene}.fasta filtered-sequences/{wildcards.gene}.fasta"

rule create_tree:
    input:
        expand("filtered-sequences/{gene}.fasta", gene=GENES)
    output:
        "RAxML_bestTree.rooted"
    shell:
        #"raxmlHPC -s {input} -m GTRCAT -n rooted -f I -t unrooted"
        "raxmlHPC -s {input} -m GTRCAT -n rooted"