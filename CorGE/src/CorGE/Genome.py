import gzip
import logging
import os
import shutil
import wget


class Genome:
    def __init__(self, name: str) -> None:
        self.name = name

    def get_name(self) -> str:
        """Returns the genome's name"""
        return self.name

    def download(self, output_fp: str):
        """Check that the file doesn't already exist then retrieves it
        Puts all files in output_fp/genomes/"""
        raise NotImplementedError

    def is_downloaded(self, genomes_fp: str, prot: bool, nucl: bool) -> bool:
        if prot and nucl:
            if os.path.exists(
                os.path.join(genomes_fp, self.name + ".faa")
            ) and os.path.exists(
                os.path.join(genomes_fp, "genomes/", self.name + ".fna")
            ):
                return True
        elif prot:
            if os.path.exists(os.path.join(genomes_fp, self.name + ".faa")):
                return True
        elif nucl:
            if os.path.exists(os.path.join(genomes_fp, self.name + ".fna")):
                return True
        else:
            return False


class AccessionGenome(Genome):
    def __init__(self, name: str, tx_id: str, url: str) -> None:
        super().__init__(name)
        self.tx_id = tx_id
        self.partial_url = url

    def download(self, output_fp: str):
        """Check that the file doesn't already exist then retrieves it
        Puts all files in output_fp/genomes/"""
        genomes_fp = os.path.join(output_fp, "genomes/")
        if not self.is_downloaded(genomes_fp, True, False):
            logging.info(f"Downloading protein-encoded genome for {self.name}")
            self.__download(genomes_fp, "_protein", ".faa")
        else:
            logging.warning(f"Found {self.name} protein-encoded genome, skipping...")
        if not self.is_downloaded(genomes_fp, False, True):
            logging.info(f"Downloading nucleotide-encoded genome for {self.name}")
            self.__download(genomes_fp, "_cds_from_genomic", ".fna")
        else:
            logging.warning(f"Found {self.name} nucleotide-encoded genome, skipping...")

    def __download(self, genomes_fp: str, name_str: str, ext: str):
        filename = f"{self.partial_url.split('/')[-1]}{name_str}{ext}.gz"
        url = f"{self.partial_url}/{filename}"

        try:
            wget.download(url, out=genomes_fp)
        except Exception as e:
            logging.error(str(e))

        with gzip.open(os.path.join(genomes_fp, filename), "rb") as f_in:
            with open(os.path.join(genomes_fp, f"{self.name}{ext}"), "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        try:
            os.remove(os.path.join(genomes_fp, filename))
        except OSError as e:
            logging.error(str(e))


class LocalGenome(Genome):
    def __init__(self, name: str, fp: str) -> None:
        super().__init__(name)
        self.fp = fp

    def download(self, output_fp: str):
        """Check that the file doesn't already exist then retrieves it
        Puts all files in output_fp/genomes/"""
        genomes_fp = os.path.join(output_fp, "genomes/")
        if not self.is_downloaded(genomes_fp, True, False):
            logging.info(f"Downloading protein-encoded genome for {self.name}")
            self.__move(genomes_fp, ".faa")
        else:
            logging.warning(f"Found {self.name} protein-encoded genome, skipping...")
        if not self.is_downloaded(genomes_fp, False, True):
            logging.info(f"Downloading nucleotide-encoded genome for {self.name}")
            self.__move(genomes_fp, ".fna")
        else:
            logging.warning(f"Found {self.name} nucleotide-encoded genome, skipping...")

    def __move(self, genomes_fp: str, ext: str):
        shutil.copyfile(
            os.path.join(self.fp, f"{self.name}{ext}"),
            os.path.join(genomes_fp, f"{self.name}{ext}"),
        )
