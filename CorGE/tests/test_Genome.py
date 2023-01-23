import pytest
from src.CorGE.Genome import Genome


@pytest.fixture
def genome():
    yield Genome("GCF_000016525.1")


def test_genome(genome):
    g = genome
    assert g.get_name() == "GCF_000016525.1"

    with pytest.raises(NotImplementedError) as e_info:
        g.download("/")

    assert g.is_downloaded("test-data/collected-genomes", True, True) == True
