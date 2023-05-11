import pytest
from src.CorGE.Genome import Genome, AccessionGenome, LocalGenome


@pytest.fixture
def genome():
    yield Genome("GCF_000016525.1")


@pytest.fixture
def accession_genome():
    yield AccessionGenome(
        "GCF_000016525.1",
        "2173",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/525/GCF_000016525.1_ASM1652v1",
    )


@pytest.fixture
def local_genome():
    yield LocalGenome(
        "GCF_000016525.1",
        "test-data/collected-genomes/",
    )


def test_genome(genome):
    g = genome
    assert g.get_name() == "GCF_000016525.1"

    with pytest.raises(NotImplementedError):
        g.download("/")


#    assert g.is_downloaded("test-data/collected-genomes", True, True) == True
#    assert g.is_downloaded("test-data/collected-genomes", True, False) == True
#    assert g.is_downloaded("test-data/collected-genomes", False, True) == True
#    assert g.is_downloaded("test-data/collected-genomes", False, False) == False


# def test_accession_genome(accession_genome):
#    g = accession_genome
#    g.download("./")
#    assert os.listdir() == ""


# def test_local_genome(local_genome):
#    g = local_genome
#    g.download("./")
#    assert os.listdir() == ""
