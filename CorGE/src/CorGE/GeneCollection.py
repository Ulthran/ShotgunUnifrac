import collections
import csv
import logging
import os
import pyhmmer.plan7
import pyhmmer.easel
import tqdm
import urllib.request
from io import TextIOWrapper


class GeneCollection:
    def __init__(
        self,
        genomes: str = os.path.join(os.getcwd(), "output/", "genomes/"),
        output: str = os.path.join(os.getcwd(), "output/"),
        file_type: str = "prot",
        name_type: str = "tx_id",
        outgroup: str = "2173",
    ) -> None:
        self.genomes = genomes
        if self.genomes[-1] != "/":
            self.genomes += "/"
        self.output = output
        if self.output[-1] != "/":
            self.output += "/"
        self.file_type = file_type
        self.name_type = name_type
        self.outgroup = outgroup

        self.filtered_fp = os.path.join(self.output, "filtered-sequences/")
        self.merged_fp = os.path.join(self.output, "merged-sequences/")
        if not os.path.exists(self.filtered_fp):
            logging.info("Making filtered sequences directory.")
            os.makedirs(self.filtered_fp)
        if not os.path.exists(self.merged_fp):
            logging.info("Making merged sequences directory.")
            os.makedirs(self.merged_fp)

        self.filtered_nucl_fp = os.path.join(self.output, "filtered-nucl-sequences/")
        if self.file_type == "nucl":
            if not os.path.exists(self.filtered_nucl_fp):
                logging.info("Making filtered nucleotide sequences directory.")
                os.makedirs(self.filtered_nucl_fp)

        self.config_fp = os.path.join(self.output, "config.yml")
        self.assembly_summary_fp = os.path.join(
            "/".join(self.genomes.split("/")[:-2]), "assembly_summary.txt"
        )

    def filter_prot(self):
        """Filters SCCGs from protein files"""
        prot_fps = [fp for fp in os.listdir(self.genomes) if fp[-4:] == ".faa"]
        prot_fps = self.__no_repeat_filter(prot_fps, self.filtered_fp)
        prot_fps = [os.path.join(self.genomes, fp) for fp in prot_fps]
        logging.debug(f"Filtering: {prot_fps}")

        logging.info("Filtering sequences...")
        with tqdm.tqdm(total=len(prot_fps)) as pbar:
            for prot_fp in prot_fps:
                pbar.update(1)
                with pyhmmer.easel.SequenceFile(prot_fp, digital=True) as seqs_file:
                    proteins = list(seqs_file)

                results = self.__run_hmmscan(proteins)

                name = prot_fp.split("/")[-1].split(".faa")[0]

                for result in results:
                    with open(
                        os.path.join(self.filtered_fp, f"{result.cog}__{name}.faa"), "w"
                    ) as f:
                        with open(prot_fp) as g:
                            self.__write_sequence(f, g, result.query)

                logging.info(f"Filtered {name}, top bitscores:")
                for result in results[:10]:
                    logging.info(
                        f"{result.query}\t{'{:.1f}'.format(result.bitscore)}\t{result.cog}"
                    )
                with open(os.path.join(self.filtered_fp, f".done_{name}"), "w") as f:
                    f.write("")

    def filter_nucl(self):
        """Filters SCCGs from nucleotide files, only runs if file_type is nucl"""
        if self.file_type == "nucl":
            nucl_fps = [
                os.path.join(self.genomes, fp)
                for fp in os.listdir(self.genomes)
                if fp[-4:] == ".fna"
            ]

            for fp in os.listdir(self.filtered_fp):
                try:
                    cog = fp.split("__")[0]
                    acc = fp.split("__")[1].split(".faa")[0]
                    matching_nucl_fp = [
                        fp for fp in nucl_fps if fp.split("/")[-1][:-4] == acc
                    ][0]

                    query = ""
                    with open(os.path.join(self.filtered_fp, fp)) as f:
                        query = f.readline().strip()[1:].split(" ")[0]

                    with open(
                        os.path.join(self.filtered_nucl_fp, f"{cog}__{acc}.fna"), "w"
                    ) as f:
                        with open(matching_nucl_fp) as g:
                            self.__write_sequence(f, g, query)
                except IndexError as e:
                    logging.debug(
                        f"File {fp} doesn't meet naming standards, skipping..."
                    )

    def merge(self):
        """Merges filtered sequences into per-SCCG files"""
        if len(os.listdir(self.merged_fp)) > 0:
            logging.warning(f"Found existing files in {self.merged_fp}, overwriting...")
            for fp in os.listdir(self.merged_fp):
                try:
                    os.remove(os.path.join(self.merged_fp, fp))
                except FileNotFoundError as e:
                    logging.debug(
                        f"While removing {os.path.join(self.merged_fp, fp)}: {e}"
                    )

        all_filtered_seqs = list()
        if self.file_type == "prot":
            all_filtered_seqs = [
                os.path.join(self.filtered_fp, fp)
                for fp in os.listdir(self.filtered_fp)
                if fp[0] != "."
            ]
        else:
            all_filtered_seqs = [
                os.path.join(self.filtered_nucl_fp, fp)
                for fp in os.listdir(self.filtered_nucl_fp)
                if fp[0] != "."
            ]

        names = self.__get_names_map(
            all_filtered_seqs, self.name_type, self.assembly_summary_fp
        )
        if self.outgroup not in names.values():
            logging.error(
                f"Didn't find outgroup {self.outgroup} in names map, make sure to use an outgroup name that matches the name you will see on the tree"
            )
            if self.outgroup in names.keys():
                logging.info(f"Did you mean {names[self.outgroup]}?")

        for fp in all_filtered_seqs:
            seq = fp.split("/")[-1]
            name = seq.split("__")[1][:-4]
            cog = seq.split("__")[0]
            new_name = names[name]
            with open(fp) as f:
                with open(os.path.join(self.merged_fp, f"{cog}.fasta"), "a") as g:
                    g.write(f"> {new_name}\n")
                    g.write(f.readlines()[1])

    def write_config(self):
        """Writes a config file for the snakemake pipeline"""
        all_merged_seqs = os.listdir(self.merged_fp)
        cogs = [fp.split(".fasta")[0] for fp in all_merged_seqs]

        cfg = "# Config file for tree building pipeline\n# Genes to build trees from\nGENES: ["
        for cog in cogs:
            cfg += f'"{cog.strip()}", '
        cfg += "]\n"

        cfg += "# Outgroup to use for rooting, false if outgroup rooting shouldn't be used\n"
        cfg += f"OUTGROUP: {self.outgroup}\n"

        cfg += f"# File type contained in merged-sequences (prot or nucl)\nTYPE: {self.file_type}\n"

        cfg += f"# Method to build tree (supermat or genetree)\nALG: supermat\n"

        cfg += f"# Directory to look for output data in\nDATA: {self.output}"

        with open(self.config_fp, "w") as f:
            f.write(cfg)

    # From pyhmmer docs https://pyhmmer.readthedocs.io/en/stable/examples/fetchmgs.html
    @staticmethod
    def __run_hmmscan(proteins: list) -> list:
        url = "https://github.com/motu-tool/fetchMGs/raw/master/lib/MG_BitScoreCutoffs.allhits.txt"

        cutoffs = {}
        with urllib.request.urlopen(url) as f:
            for line in csv.reader(TextIOWrapper(f), dialect="excel-tab"):
                if not line[0].startswith("#"):
                    cutoffs[line[0]] = float(line[1])

        baseurl = "https://github.com/motu-tool/fetchMGs/raw/master/lib/{}.hmm"

        hmms = []
        for cog in cutoffs:
            with urllib.request.urlopen(baseurl.format(cog)) as f:
                hmm = next(pyhmmer.plan7.HMMFile(f))
                cutoff = cutoffs[hmm.name.decode()]
                hmm.cutoffs.trusted = (cutoff, cutoff)
                hmms.append(hmm)

        Result = collections.namedtuple("Result", ["query", "cog", "bitscore"])

        results = []
        for top_hits in pyhmmer.hmmsearch(hmms, proteins, bit_cutoffs="trusted"):
            for hit in top_hits:
                cog = hit.best_domain.alignment.hmm_name.decode()
                results.append(Result(hit.name.decode(), cog, hit.score))

        best_results = {}
        keep_query = set()
        for result in results:
            if result.query in best_results:
                previous_bitscore = best_results[result.query].bitscore
                if result.bitscore > previous_bitscore:
                    best_results[result.query] = result
                    keep_query.add(result.query)
                elif result.bitscore == previous_bitscore:
                    if best_results[result.query].cog != hit.cog:
                        keep_query.remove(result.query)
            else:
                best_results[result.query] = result
                keep_query.add(result.query)

        return [best_results[k] for k in sorted(best_results) if k in keep_query]

    @staticmethod
    def __no_repeat_filter(prot_fps: list, filtered_fp: str) -> list:
        filtered_prot_fps = []
        for fp in prot_fps:
            existing_filtered_files = [
                filtered_fp
                for filtered_fp in os.listdir(filtered_fp)
                if fp in filtered_fp
            ]
            l = len(existing_filtered_files)
            if os.path.exists(os.path.join(filtered_fp, f".done_{fp[:-4]}")):
                logging.info(
                    f"Skipping protein filter step for {fp[:-4]} because all files already exist..."
                )
            elif l > 0 and not os.path.exists(
                os.path.join(filtered_fp, f".done_fp[:-4]")
            ):
                logging.warning(
                    f"Found partial protein filter files for {fp[:-4]}, overwriting..."
                )
                filtered_prot_fps.append(fp)
            else:
                filtered_prot_fps.append(fp)

        return filtered_prot_fps

    @staticmethod
    def __write_sequence(
        out: TextIOWrapper, seqs: TextIOWrapper, query: str, acc: str = None
    ):
        seq = ""
        add = False
        for l in seqs.readlines():
            if add:
                if l[0] != ">":
                    seq += l.strip()
                else:
                    break
            if query in l and l[0] == ">":
                add = True
                if acc:
                    seq += f"> {acc}\n"
                else:
                    seq += f"{l.strip()}\n"

        out.write(seq + "\n")

    @staticmethod
    def __get_names_map(fps: list, nt: str, assembly_summary_fp: str) -> dict:
        accs = {
            fp.split("/")[-1]
            .split("__")[1]
            .split(".f")[0]: fp.split("/")[-1]
            .split("__")[1]
            .split(".f")[0]
            for fp in fps
        }
        header = ""
        if nt == "acc":
            logging.debug(f"Names map: {accs}")
            return accs
        elif nt == "tx_id":
            header = "species_taxid"
        elif nt == "strain":
            header = "infraspecific_name"
        elif nt == "species":
            header = "organism_name"
        accs_list = list(accs.keys())
        ids = {}

        with open(assembly_summary_fp) as f:
            reader = csv.reader(f, dialect=csv.excel_tab)
            next(reader)  # First row is a comment
            headers = next(reader)  # This row has the headers
            headers[0] = headers[0][
                2:
            ]  # Remove the "# " from the beginning of the first element

            idIndex = headers.index(header)
            accIndex = headers.index("assembly_accession")

            for line in reader:
                try:
                    if line[accIndex] in accs_list:
                        ids[line[accIndex]] = line[idIndex]
                except IndexError:
                    None  # Incomplete assembly_summary entry

        # Fill in names that didn't find a match
        for name in [name for name in accs.keys() if name not in ids.keys()]:
            ids[name] = name

        logging.debug(f"Names map: {ids}")
        return ids
