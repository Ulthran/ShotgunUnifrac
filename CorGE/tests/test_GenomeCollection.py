import os
import pytest
from src.CorGE.GenomeCollection import GenomeCollection

OUTPUT_FP = os.path.join(os.getcwd(), "output/")


@pytest.fixture
def genome_collection():
    yield GenomeCollection()


def test_genome_collection(genome_collection):
    gc = genome_collection
    gc.dryrun()
    gc.collect()


@pytest.fixture
def genome_collection_species():
    yield GenomeCollection(OUTPUT_FP, [2173])


def test_genome_collection(genome_collection):
    gc = genome_collection
    gc.dryrun()
    gc.collect()


@pytest.fixture
def genome_collection_accessions():
    yield GenomeCollection(OUTPUT_FP, ["GCF_000007725.1"])


def test_genome_collection(genome_collection):
    gc = genome_collection
    gc.dryrun()
    gc.collect()


@pytest.fixture
def genome_collection_local():
    yield GenomeCollection(OUTPUT_FP, "test-data/collected-genomes")


def test_genome_collection(genome_collection):
    gc = genome_collection
    gc.dryrun()
    gc.collect()
