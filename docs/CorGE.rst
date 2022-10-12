.. _CorGE:

=====================
CorGE Usage Guide
=====================

.. contents::
   :depth: 2

Options
*******

CorGE
-----

.. code-block:: shell

    CorGE [-h] {collect_genomes,extract_genes} ...

.. code-block:: shell

    -h/--help: Display help.

CorGE collect_genomes
---------------------

Collect nucleotide- and protein-encoded genomes of interest.

.. code-block:: shell

    CorGE collect_genomes [-h] [--all] [--ncbi_species NCBI_SPECIES] [--ncbi_accessions NCBI_ACCESSIONS] [--local LOCAL] [--outgroup OUTGROUP] [-n] output_dir

.. code-block:: shell

    -h/--help: Display help.
    --all: Collect one representative genome from each species listed in NCBI's RefSeq database. Don't use this with --ncbi_species, --ncbi_accessions, or --local.
    --ncbi_species: File listing species level taxon ids to be collected from NCBI.
    --ncbi_accessions: File listing genome accessions to be collected from NCBI.
    --local: Directory containing nucleotide- and protein-encoded pairs of genome files. Any unpaired files will be ignored.
    --outgroup: Specify the outgroup for tree rooting. Integers will be parsed as species level taxon ids and retrieved from NCBI. Otherwise will search for a matching nucleotide-encoded file in ouput_dir or local (Default: 2173, enter None to not use outgroup rooting).
    -n: Dry run, show what would be gathered but don't do it.
    output_dir: Directory to collect genomes in.

CorGE extract_genes
-------------------

Extract SCCGs from all collected genomes and curate data for tree building.

.. code-block:: shell

    CorGE extract_genes [-h] [-o OUTPUT] [-t {prot,nucl}] [-n {acc,txid,strain,species}] genomes

.. code-block:: shell

    -h/--help: Display help.
    -o/--output: Directory to write output to (Default: ./).
    -t/--type: Output in merged-sequences can be nucleotide- or protein-encoded (Default: prot).
    -n/--name: Names to show on final tree (Default: txid).
    genomes: Directory with collected genomes (curated with collect_genomes).

collect_genomes
***************

This is the command for retrieving genomes for all of the bacteria you want in your tree, from NCBI, a local directory, or both. `CorGE` does all of its work with pairs of files: a protein-encoded genome (usually saved with a `.faa` extension) and a nucleotide-encoded genome (usually saved with a `.fna` extension). These files have the same name but different extensions with the name for any NCBI files being the genome accession (e.g. "GCF_000218545") and the name for any local files being preserved.

.. tip::

    Especially when downloading lots of genomes, it's always a good idea to run your command with the dryrun option (`-n`) first. This will print out what it's planning to do without doing it, so you can verify that it will do the right thing.

Common Use Cases
----------------

To collect a set of representative/reference genomes for a list of species, create a file of species-level taxon ids, one id per line, and pass that file to `CorGE`:

.. code-block:: shell

    CorGE collect_genomes . --ncbi_species species_list.txt

.. tip::

    If you don't care about rooting the final tree, you can specify `--outgroup None`. The pipeline will still use a midpoint algorithm to root the final tree, but the input to that step will be the unrooted tree. 

To collect one genome for each species NCBI has, use the `--all` option:

.. code-block:: shell

    CorGE collect_genomes /path/to/db --all

Suppose you want to create a strain level tree from some existing NCBI E. coli genomes and some that you have locally and then root that tree with a reference Clostridium botulinum genome. You create a list of the genome accessions you want to collect (same as species taxa, one accession per line in a text file) and run:

.. code-block:: shell

    CorGE collect_genomes ecoli-db/ --ncbi_accessions accession_list.txt

This will collect each of those genomes and put them in `ecoli-db/` (as well as grabbing Methanobrevibacter smithii in the `outgroup` dir, this will be overwritten by the next step). Next we need to get all the local files in there, but we need them to follow a couple rules: 1) only pairs of files will be collected so every nucleotide file should have a paired protein file with the same name (e.g. `example.fna` and `example.faa`) and 2) the annotations for the protein files should have a unique name as the first space-separated piece of the annotation and then corresponding sequences in the nucleotide files should have annotations that contain that name (e.g. first protein sequence is annotated with `> example00001 description and so on` and corresponding nucleotide sequence is annotated with `> WXK40_example00001_DEADBEEF`). Provided they are in compliance with the above and all in one directory run:

.. code-block:: shell

    CorGE collect_genomes ecoli-db/ --local assembled-genomes/ --outgroup 1491

This should leave you with all your files organized into `ecoli-db/protein`, `ecoli-db/nucleotide`, and `ecoli-db/outgroup` directories.

.. tip::

    You could also do this all in one command `CorGE collect_genomes ecoli-db/ --ncbi_accessions accession_list.txt --local assembled-genomes/ --outgroup 1491`

extract_genes
*************

To curate genes for a multi-species bacteria tree run:

.. code-block:: shell

    CorGE extract_genes /path/to/db

Continuing the example from the last section, to curate genes for this strain level E. coli tree, run:

.. code-block:: shell

    CorGE extract_genes ecoli-db/ --type nucl --name strain

This will prepare you to build a nucleotide based tree of E. coli strains where the leaf names will be either the strain name (if it came from NCBI) or the file name (if it was local).