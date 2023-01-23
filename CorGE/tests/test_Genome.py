import pytest
from src.CorGE.Genome import Genome, AccessionGenome, LocalGenome


@pytest.fixture
def genome():
    yield Genome("GCF_000016525.1")

@pytest.fixture
def accession_genome(name, txid, url):
    yield AccessionGenome(name, txid, url)


def test_genome(genome):
    g = genome
    assert g.get_name() == "GCF_000016525.1"

    with pytest.raises(NotImplementedError) as e_info:
        g.download("/")

    assert g.is_downloaded("test-data/collected-genomes", True, True) == True
    assert g.is_downloaded("test-data/collected-genomes", True, False) == True
    assert g.is_downloaded("test-data/collected-genomes", False, True) == True
    assert g.is_downloaded("test-data/collected-genomes", False, False) == False

def test_accession_genome_txid(accession_genome):
    g = accession_genome("GCF_000016525.1", "2173")