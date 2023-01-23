import pytest
from src.CorGE.Genome import Genome


@pytest.fixture
def genome():
    yield Genome("TEST")


def test_genome(genome):
    g = genome
    assert g.get_name() == "TEST"
    with pytest.raises(NotImplementedError) as e_info:
        g.download("/")
