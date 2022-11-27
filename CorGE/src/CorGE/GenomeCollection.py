import logging
import os
from .Genome import SpeciesGenome, AccessionGenome, LocalGenome

class GenomeCollection():
    def __init__(self,
                    output_fp: str = os.path.join(os.getcwd(), "projects/"),
                    species: list = list(),
                    accessions: list = list(),
                    local_fp: str = None,
                    outgroup: str = "2173") -> None:
        self.genomes = list()
        for s in species:
            self.genomes.append(SpeciesGenome(s))
        for a in accessions:
            self.genomes.append(AccessionGenome(a))
        for l in self.list_valid_genomes(local_fp):
            self.genomes.append(LocalGenome(l, local_fp))

        self.output_fp = output_fp
        if not os.path.exists(self.output_fp):
            logging.info("Making output directory.")
            os.makedirs(self.output_fp)

        self.filtered_fp = os.path.join(self.output_fp, "filtered-sequences/")
        self.merged_fp = os.path.join(self.output_fp, "merged-sequences/")
        if not os.path.exists(self.filtered_fp):
            logging.info("Making filtered sequences directory.")
            os.makedirs(self.filtered_fp)
        if not os.path.exists(self.merged_fp):
            logging.info("Making merged sequences directory.")
            os.makedirs(self.merged_fp)
        
        self.outgroup = outgroup
        if self.outgroup.isdigit():
            logging.info("Interpreting outgroup as species-level taxon id.")
            self.genomes.append(SpeciesGenome(self.outgroup))
        elif os.path.exists(os.path.join(local_fp, self.outgroup)):
            logging.info(f"Interpreting outgroup as name in {local_fp}.")
            self.genomes.append(LocalGenome(self.outgroup, local_fp))
        else:
            logging.info("Interpreting outgroup as an NCBI accession.")
            self.genomes.append(AccessionGenome(self.outgroup))

    def dryrun(self):
        """Prints all the accessions and filepaths to be collected"""
        accs = [g.get_accession_or_fp() for g in self.genomes if type(g) == AccessionGenome or type(g) == SpeciesGenome]
        fps = [g.get_accession_or_fp() for g in self.genomes if type(g) == AccessionGenome]
        logging.info(f"Accessions to be collected (total: {len(accs)}): {accs}")
        logging.info(f"Local genomes to be collected (total: {len(fps)}): {fps}")
    
    def collect(self):
        """Downloads or copies all Genome objects in genomes to output_fp"""
        for g in self.genomes:
            g.download(self.output_fp)

    def filter_prot(self):
        """Filters SCCGs from protein files"""
        prot_fps = [os.path.join(self.output_fp, fp) for fp in os.listdir(self.output_fp) if fp[-4:] == ".faa"]

    def filter_nucl(self):
        """Filters SCCGs from nucleotide files"""

    def merge(self):
        """Merges filtered sequences into per-SCCG files"""

    def write_config(self):
        """Writes a config file for the snakemake pipeline"""

    def __update_collection(self):
        """Updates genomes list to match the current output_fp/genomes/ dir
        This should be used to allow for easy continuation curating an existing data set"""

    @staticmethod
    def list_valid_genomes(fp: str) -> list:
        """Returns a list of filenames that have both .faa and .fna files in the given filepath"""
        if not os.path.exists(fp):
            logging.warning(f"Path {fp} doesn't exist, skipping local genomes.")
            return list()
        files = os.listdir(fp)
        return [f[-4:] for f in files if f[-4:] == '.fna' and f"{f[-4:]}.faa" in files]

    @staticmethod
    def is_protein_file(fp: str) -> bool:
        """Tells whether or not an input file is a protein file
        (relative path is assumed to be from output_fp/genomes/)"""
        #TODO: Account for gzip files
        if os.path.split(fp)[1][-4:] == ".faa":
            return True
        return False