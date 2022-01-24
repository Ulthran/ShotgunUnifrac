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
IDS_FULL = [
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
IDS = [
    "GCF_000196555.1_ASM19655v1",
    "GCF_902387545.1_UHGG_MGYG-HGUT-02478",
    "GCF_020885855.1_ASM2088585v1",
    "GCF_020735445.1_ASM2073544v1"
]
GENES = ["secE", "secG"]

rule all:
    input: "rooted.tree"

rule all_gene_files:
    input: expand("sequences/{gene}.fasta", gene=GENES)

rule get_ncbi_sequences:
    input:
        "data/sample.fastq"
    output:
        expand("ncbi/{id}_cds_from_genomic.fasta", id=IDS)
    shell:
        "./download-genes.sh {output}"

rule extract_marker_genes:
    input:
        "ncbi/{id}_cds_from_genomic.fasta"
    output:
        expand("sequences/{gene}__{{id}}.fasta", gene=GENES)
    shell:
        "python3 -c 'import filterGenes; filterGenes.filter_seq_genes_sm(\"ncbi/{wildcards.id}_cds_from_genomic.fasta\", {output})' > {output}"

rule merge_marker_genes:
    input:
        expand("sequences/{{gene}}__{id}.fasta", id=IDS)
    output:
        temp("sequences/{gene}.fasta")
    shell:
        "cat sequences/{wildcards.gene}__*.fasta > {output}"

rule align_fasta:
    input:
        "sequences/{gene}.fasta"
    output:
        temp("aligned-sequences/{gene}.fasta")
    shell:
        "muscle sequences/{gene}.fasta --output aligned-sequences/{gene}.fasta"

rule filter_columns:
    # Identity op for now, okfasta
    input:
        "aligned-sequences/{gene}.fasta"
    output:
        temp("filtered-sequences/{gene}.fasta")
    shell:
        "cp aligned-sequences/{gene}.fasta filtered-sequences/{gene}.fasta"

rule create_tree:
    input:
        expand("filtered-sequences/{gene}.fasta", gene=GENES)
    output:
        "rooted.tree"
    shell:
        "raxml filtered-sequences/{gene}.fasta -f I"