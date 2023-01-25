import os
import pytest
import shutil
from src.CorGE.GeneCollection import GeneCollection
from . import *


@pytest.fixture
def gene_collection():
    yield GeneCollection(
        os.path.join(TEST_DATA_FP, "collected-genomes"), OUTPUT_FP, "nucl"
    )


def test_gene_collection(gene_collection):
    gc = gene_collection
    gc.filter_prot()
    assert set(
        [
            "COG0012__GCF_000016525.1.faa",
            "COG0012__GCF_000007725.1.faa",
            "COG0541__GCF_000016525.1.faa",
        ]
    ).issubset(set(os.listdir(FILTERED_FP)))
    gc.filter_nucl()
    assert set(
        [
            "COG0012__GCF_000016525.1.fna",
            "COG0012__GCF_000007725.1.fna",
            "COG0541__GCF_000016525.1.fna",
        ]
    ).issubset(set(os.listdir(FILTERED_NUCL_FP)))
    gc.merge()
    assert len(os.listdir(MERGED_FP)) == 41
    gc.write_config()
    assert "config.yml" in os.listdir(OUTPUT_FP)
