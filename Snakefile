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
GENES = ["secE", "secG"]

rule all:
    input: "trees/RAxML_rootedTree.final"

rule all_gene_files:
    input: expand("sequences/{gene}.fasta", gene=GENES)

rule all_filtered_files:
    input: expand("filtered-sequences/{gene}.fasta", gene=GENES)

rule get_ncbi_sequences:
    input:
        "data/sample.fastq"
    output:
        temp("ncbi/{id}_cds_from_genomic.fasta")
    log:
        "logs/get_ncbi_sequences/{id}.log"
    shell:
        "./download-genes.sh {output}"

rule extract_marker_genes:
    input:
        "ncbi/{id}_cds_from_genomic.fasta"
    output:
        temp("sequences/{gene}__{id}.fasta")
    log:
        "logs/extract_marker_genes/{gene}__{id}.log"
    shell:
        "python3 -c 'import filterGenes; filterGenes.filter_seq_genes_sm(\"ncbi/{wildcards.id}_cds_from_genomic.fasta\", \"{output}\")' > {output}"

rule merge_marker_genes:
    input:
        expand("sequences/{gene}__{id}.fasta", id=IDS, gene=GENES)
    output:
        "sequences/{gene}.fasta"
    log:
        "logs/merge_marker_genes/{gene}.log"
    shell:
        # Rename sequences to human-readable leaf names
        "cat sequences/{wildcards.gene}__*.fasta > {output}"

rule align_fasta:
    input:
        "sequences/{gene}.fasta"
    output:
        temp("aligned-sequences/{gene}.fasta")
    log:
        "logs/align_fasta/{gene}.log"
    shell:
        "muscle -in sequences/{wildcards.gene}.fasta -out aligned-sequences/{wildcards.gene}.fasta"

rule filter_columns:
    # Identity op for now, okfasta
    input:
        "aligned-sequences/{gene}.fasta"
    output:
        "filtered-sequences/{gene}.fasta"
    log:
        "logs/filter_columns/{gene}.log"
    shell:
        "cp aligned-sequences/{wildcards.gene}.fasta filtered-sequences/{wildcards.gene}.fasta"

rule create_trees:
    input:
        "filtered-sequences/{gene}.fasta"
    output:
        "trees/RAxML_bestTree.{gene}"
    log:
        "logs/create_trees/RAxML_bestTree.{gene}.log"
    shell:
        "cd trees/ && "
        "raxmlHPC -s ../{input} -m GTRCAT -n {wildcards.gene} && "
        "cd .."

rule merge_trees:
    input:
        expand("trees/RAxML_bestTree.{gene}", gene=GENES)
    output:
        "trees/RAxML_rootedTree.final"
    log:
        "logs/merge_trees/final.rooted.log"
    shell:
        "cat {input} > trees/merged.in && "
        "java -jar Astral/astral.5.7.8.jar -i trees/merged.in -o trees/final.unrooted 2>logs/merge_trees/final.unrooted.out.log && "
        "cd trees/ && "
        "raxmlHPC -f I -m GTRCAT -t final.unrooted -n final && "
        "cd .."