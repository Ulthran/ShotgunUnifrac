import pytest
from src.CorGE.GenomeCollection import GenomeCollection


@pytest.fixture
def genome_collection():
    yield GenomeCollection()

def test_genome_collection(genome_collection):
    gc = genome_collection
    gc.dryrun()
    gc.collect()

