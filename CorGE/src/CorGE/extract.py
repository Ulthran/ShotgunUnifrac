from .GeneCollection import GeneCollection


def extract_genes(args: dict):
    if args["file_type"]:
        args["file_type"] = str(args["file_type"])
    if args["name_type"]:
        args["name_type"] = str(args["name_type"])

    gc_args = {
        k: v
        for k, v in args.items()
        if v and k in ["genomes", "output", "file_type", "name_type", "outgroup"]
    }

    gc = GeneCollection(**gc_args)

    gc.filter_prot()
    gc.filter_nucl()
    gc.merge()
    gc.write_config()
