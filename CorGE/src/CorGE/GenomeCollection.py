import csv
import logging
import os
import wget
from .Genome import AccessionGenome, LocalGenome


class GenomeCollection:
    def __init__(
        self,
        output_fp: str = os.path.join(os.getcwd(), "output/"),
        ncbi_species: list = [],
        ncbi_accessions: list = [],
        local_fp: str = "",
    ) -> None:
        self.output_fp = output_fp
        if self.output_fp[-1] != "/":
            self.output_fp += "/"
        if not os.path.exists(self.output_fp):
            logging.info("Making output directory.")
            os.makedirs(self.output_fp)
        logging.info(f"Output location: {self.output_fp}")
        self.assembly_summary_fp = os.path.join(self.output_fp, "assembly_summary.txt")
        if not os.path.exists(self.assembly_summary_fp):
            self.__download_assembly_summary()

        self.local_fp = (
            local_fp if os.path.isabs(local_fp) else os.path.join(os.getcwd(), local_fp)
        )

        self.genomes = list()
        for vals in self.__accessions_for_species(
            ncbi_species, self.assembly_summary_fp
        ):
            self.genomes.append(AccessionGenome(*vals))
        for vals in self.__species_for_accessions(
            ncbi_accessions, self.assembly_summary_fp
        ):
            self.genomes.append(AccessionGenome(*vals))
        for l in self.__list_valid_genomes(self.local_fp):
            self.genomes.append(LocalGenome(l, self.local_fp))

        self.genomes_fp = os.path.join(self.output_fp, "genomes/")
        if not os.path.exists(self.genomes_fp):
            logging.info("Making genomes directory.")
            os.makedirs(self.genomes_fp)

    def dryrun(self):
        """Prints all the accessions and filepaths to be collected"""
        accs = [g.get_name() for g in self.genomes if isinstance(g, AccessionGenome)]
        fps = [g.get_name() for g in self.genomes if isinstance(g, LocalGenome)]
        accs_len = len(accs)
        fps_len = len(fps)
        if accs_len > 50:
            logging.debug(f"Full accessions output: {accs}")
            accs = accs[:50]
        if fps_len > 50:
            logging.debug(f"Full filepaths output: {fps}")
            fps = fps[:50]
        logging.info(
            f"Accessions to be collected (total: {accs_len}): {accs + ['...'] if accs_len > 50 else accs}"
        )
        logging.info(
            f"Local genomes to be collected (total: {fps_len}): {fps + ['...'] if fps_len > 50 else fps}"
        )

    def all_species(self):
        """Adds one representative of every species to genomes list"""
        ids = []
        with open(self.assembly_summary_fp) as f:
            reader = csv.reader(f, dialect=csv.excel_tab)
            next(reader)  # First row is a comment
            headers = next(reader)  # This row has the headers
            headers[0] = headers[0][
                2:
            ]  # Remove the "# " from the beginning of the first element

            idIndex = headers.index("species_taxid")
            accIndex = headers.index("assembly_accession")
            refSeqIndex = headers.index("refseq_category")
            urlIndex = headers.index("ftp_path")

            for line in reader:
                try:
                    if int(line[idIndex]) not in ids and line[refSeqIndex] != "na":
                        ids.append(int(line[idIndex]))
                        self.genomes.append(
                            AccessionGenome(
                                line[accIndex], line[idIndex], line[urlIndex]
                            )
                        )
                except IndexError:
                    pass  # Incomplete assembly_summary entry

    def collect(self):
        """Downloads or copies all Genome objects in genomes to output_fp"""
        for g in self.genomes:
            g.download(self.output_fp)

    ### Private Methods

    def __download_assembly_summary(self):
        """Downloads assembly_summary.txt to output_fp"""
        logging.info("assembly_summary.txt not found, fetching...")
        wget.download(
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt",
            out=self.output_fp,
        )

        with open(self.assembly_summary_fp, "a") as f:  # Append 2173 default outgroup
            f.write(
                "GCF_000016525.1\tPRJNA224116\tSAMN02604313\t\trepresentative genome\t420247\t2173\tMethanobrevibacter smithii ATCC 35061\tstrain=ATCC 35061; PS; DSMZ 861\t\tlatest\tComplete Genome\tMajor\tFull\t2007/06/04\tASM1652v1\tWashington University Center for Genome Sciences\tGCA_000016525.1\tidentical\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/525/GCF_000016525.1_ASM1652v1\t\tassembly from type material\tna"
            )

    @staticmethod
    def __accessions_for_species(species: list, assembly_summary: str) -> list:
        """Return list of one accession for each taxon id"""
        ret = list()
        tmp_species = species
        with open(assembly_summary) as f:
            reader = csv.reader(f, dialect=csv.excel_tab)
            next(reader)  # First row is a comment
            headers = next(reader)  # This row has the headers
            headers[0] = headers[0][
                2:
            ]  # Remove the "# " from the beginning of the first element

            idIndex = headers.index("species_taxid")
            accIndex = headers.index("assembly_accession")
            refSeqIndex = headers.index("refseq_category")
            urlIndex = headers.index("ftp_path")

            for line in reader:
                try:
                    if line[idIndex] in tmp_species and line[refSeqIndex] != "na":
                        ret.append((line[accIndex], line[idIndex], line[urlIndex]))
                        tmp_species.remove(line[idIndex])
                except IndexError:
                    pass  # Incomplete assembly_summary entry

        return ret

    @staticmethod
    def __species_for_accessions(accs: list, assembly_summary: str) -> list:
        """Return list of one taxon id for each accession"""
        ret = list()
        with open(assembly_summary) as f:
            reader = csv.reader(f, dialect=csv.excel_tab)
            next(reader)  # First row is a comment
            headers = next(reader)  # This row has the headers
            headers[0] = headers[0][
                2:
            ]  # Remove the "# " from the beginning of the first element

            idIndex = headers.index("species_taxid")
            accIndex = headers.index("assembly_accession")
            urlIndex = headers.index("ftp_path")

            for line in reader:
                try:
                    if line[accIndex] in accs:
                        ret.append((line[accIndex], line[idIndex], line[urlIndex]))
                except IndexError:
                    None  # Incomplete assembly_summary entry

        return ret

    @staticmethod
    def __list_valid_genomes(fp: str) -> list:
        """Returns a list of filenames that have both .faa and .fna files in the given filepath"""
        if fp:
            if not os.path.exists(fp):
                logging.warning(f"Path {fp} doesn't exist, skipping local genomes.")
                return []
            files = os.listdir(fp)
            return [
                f[:-4] for f in files if f[-4:] == ".fna" and f"{f[:-4]}.faa" in files
            ]
        return []

    @staticmethod
    def __is_protein_file(fp: str) -> bool:
        """Tells whether or not an input file is a protein file
        (relative path is assumed to be from output_fp/genomes/)"""
        # TODO: Account for gzip files
        if os.path.split(fp)[1][-4:] == ".faa":
            return True
        return False
