import os
import pytest
import shutil
from src.CorGE.GenomeCollection import GenomeCollection
from . import *


@pytest.fixture
def genome_collection():
    yield GenomeCollection(OUTPUT_FP)
    shutil.rmtree(OUTPUT_FP)


def test_genome_collection(genome_collection):
    gc = genome_collection
    gc.dryrun()
    gc.collect()

    assert set(os.listdir(GENOMES_FP)) == set([])


@pytest.fixture
def genome_collection_species():
    yield GenomeCollection(OUTPUT_FP, ["2173"])
    shutil.rmtree(OUTPUT_FP)


def test_genome_collection_species(genome_collection_species):
    gc = genome_collection_species
    gc.dryrun()
    gc.collect()

    assert set(os.listdir(GENOMES_FP)) == set(
        ["GCF_000016525.1.fna", "GCF_000016525.1.faa"]
    )


@pytest.fixture
def genome_collection_accessions():
    yield GenomeCollection(OUTPUT_FP, [], ["GCF_000007725.1"])
    shutil.rmtree(OUTPUT_FP)


def test_genome_collection_accessions(genome_collection_accessions):
    gc = genome_collection_accessions
    gc.dryrun()
    gc.collect()

    assert set(os.listdir(GENOMES_FP)) == set(
        ["GCF_000007725.1.fna", "GCF_000007725.1.faa"]
    )


@pytest.fixture
def genome_collection_local():
    yield GenomeCollection(
        OUTPUT_FP, [], [], os.path.join(TEST_DATA_FP, "collected-genomes")
    )
    shutil.rmtree(OUTPUT_FP)


def test_genome_collection_local(genome_collection_local):
    gc = genome_collection_local
    gc.dryrun()
    gc.collect()

    assert set(os.listdir(GENOMES_FP)) == set(
        [
            "GCF_000007725.1.faa",
            "GCF_000007725.1.fna",
            "GCF_000010525.1.faa",
            "GCF_000010525.1.fna",
            "GCF_000012885.1.faa",
            "GCF_000012885.1.fna",
            "GCF_000016525.1.faa",
            "GCF_000016525.1.fna",
            "GCF_000020965.1.faa",
            "GCF_000020965.1.fna",
            "GCF_000218545.1.faa",
            "GCF_000218545.1.fna",
            "GCF_000378225.1.faa",
            "GCF_000378225.1.fna",
            "GCF_001375595.1.faa",
            "GCF_001375595.1.fna",
            "GCF_001735525.1.faa",
            "GCF_001735525.1.fna",
            "GCF_007197645.1.faa",
            "GCF_007197645.1.fna",
            "GCF_023159115.1.faa",
            "GCF_023159115.1.fna",
            "GCF_900111765.1.faa",
            "GCF_900111765.1.fna",
        ]
    )
