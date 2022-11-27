import os

class Genome():
    def __init__(self, name: str) -> None:
        self.name = name
    
    def get_name(self) -> str:
        """Returns the genome's name"""
        return self.name

    def get_accession_or_fp(self) -> str:
        """Returns the genome's accession or filepath (for local genomes)"""
        raise NotImplementedError

    def url_for(self) -> str:
        """Returns the NCBI URL for the genome (without filename)"""
        raise NotImplementedError("")
    
    def download(self, output_fp: str, nucleotide: bool = True, protein: bool = True):
        """Check that the file doesn't already exist then retrieves it
        Puts all files in output_fp/genomes/"""
        raise NotImplementedError("")
    
    def get_protein_fp(self) -> str:
        """Returns the filepath to the protein-encoded genome or None if not downloaded"""
        raise NotImplementedError("")
    
    def get_nucleotide_fp(self):
        """Returns the filepath to the nucleotide-encoded genome or None if not downloaded"""
        raise NotImplementedError("")

class SpeciesGenome(Genome):
    def __init__(self, name: str) -> None:
        super().__init__(name)

    def get_accession_or_fp(self) -> str:
        """Returns the genome's accession or filepath (for local genomes)"""
        raise NotImplementedError

    def url_for(self) -> str:
        """Returns the NCBI URL for the genome (without filename)"""
        raise NotImplementedError("")
    
    def download(self, output_fp: str, nucleotide: bool = True, protein: bool = True):
        """Check that the file hasn't already been downloaded then retrieves it"""
        raise NotImplementedError("")
    
    def get_protein_fp(self) -> str:
        """Returns the filepath to the protein-encoded genome or None if not downloaded"""
        raise NotImplementedError("")
    
    def get_nucleotide_fp(self):
        """Returns the filepath to the nucleotide-encoded genome or None if not downloaded"""
        raise NotImplementedError("")

class AccessionGenome(Genome):
    def __init__(self, name: str) -> None:
        super().__init__(name)

    def get_accession_or_fp(self) -> str:
        """Returns the genome's accession or filepath (for local genomes)"""
        raise NotImplementedError

    def url_for(self) -> str:
        """Returns the NCBI URL for the genome (without filename)"""
        raise NotImplementedError("")
    
    def download(self, output_fp: str, nucleotide: bool = True, protein: bool = True):
        """Check that the file hasn't already been downloaded then retrieves it"""
        raise NotImplementedError("")
    
    def get_protein_fp(self) -> str:
        """Returns the filepath to the protein-encoded genome or None if not downloaded"""
        raise NotImplementedError("")
    
    def get_nucleotide_fp(self):
        """Returns the filepath to the nucleotide-encoded genome or None if not downloaded"""
        raise NotImplementedError("")

class LocalGenome(Genome):
    def __init__(self, name: str, fp: str) -> None:
        super().__init__(name)

    def get_accession_or_fp(self) -> str:
        """Returns the genome's accession or filepath (for local genomes)"""
        raise NotImplementedError

    def url_for(self) -> str:
        """Returns the NCBI URL for the genome (without filename)"""
        raise NotImplementedError("")
    
    def download(self, output_fp: str, nucleotide: bool = True, protein: bool = True):
        """Check that the file hasn't already been downloaded then retrieves it"""
        raise NotImplementedError("")
    
    def get_protein_fp(self) -> str:
        """Returns the filepath to the protein-encoded genome or None if not downloaded"""
        raise NotImplementedError("")
    
    def get_nucleotide_fp(self):
        """Returns the filepath to the nucleotide-encoded genome or None if not downloaded"""
        raise NotImplementedError("")

