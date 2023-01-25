from .GenomeCollection import GenomeCollection


def collect_genomes(args: dict):
    try:
        args["ncbi_species"] = [
            line.rstrip() for line in args["ncbi_species"].readlines()
        ]
    except AttributeError as e:
        None
    try:
        args["ncbi_accessions"] = [
            line.rstrip() for line in args["ncbi_accessions"].readlines()
        ]
    except AttributeError as e:
        None

    gc_args = {
        k: v
        for k, v in args.items()
        if v and k in ["output_fp", "ncbi_species", "ncbi_accessions", "local_fp"]
    }

    gc = GenomeCollection(**gc_args)

    if args["all"]:
        gc.all_species()
    if args["n"]:
        gc.dryrun()
    else:
        gc.collect()
