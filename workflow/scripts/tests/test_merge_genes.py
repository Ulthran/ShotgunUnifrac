import os
import shutil
import tempfile
import unittest

from collections import Counter
from genetools.mergeGenes import merge_genes

class MergeGenesTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        os.mkdir(os.path.join(self.temp_dir, "output"))
        os.mkdir(os.path.join(self.temp_dir, "output/sequences"))
        os.mkdir(os.path.join(self.temp_dir, "output/merged-sequences"))

        self.gene_counter = Counter({'TEST_GENE': 4, 'FAILED_TEST_GENE': 3})

        self.test_gene_1_fp = os.path.join(self.temp_dir, "output/sequences/TEST_GENE__TEST_1.fasta")
        self.test_gene_1_contents = (
            "TEST_1\n"
            "ATCG\n")
        with open(self.test_gene_1_fp, "w") as f:
            f.write(self.test_gene_1_contents)

        self.test_gene_2_fp = os.path.join(self.temp_dir, "output/sequences/TEST_GENE__TEST_2.fasta")
        self.test_gene_2_contents = (
            "TEST_2\n"
            "GCTA\n")
        with open(self.test_gene_2_fp, "w") as f:
            f.write(self.test_gene_2_contents)

        self.test_gene_3_fp = os.path.join(self.temp_dir, "output/sequences/FAILED_TEST_GENE__TEST_1.fasta")
        self.test_gene_3_contents = (
            "TEST_3\n"
            "SHOULD FAIL\n")
        with open(self.test_gene_3_fp, "w") as f:
            f.write(self.test_gene_3_contents)
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_merge_genes(self):
        merge_genes(self.gene_counter, self.temp_dir)

        with open(os.path.join(self.temp_dir, "output/merged-sequences/TEST_GENE.fasta")) as f:
            self.assertEqual(next(f), "TEST_2\n")
            self.assertEqual(next(f), "GCTA\n")
            self.assertEqual(next(f), "TEST_1\n")
            self.assertEqual(next(f), "ATCG\n")
        
        self.assertEqual(os.path.exists(os.path.join(self.temp_dir, "output/merged-sequences/FAILED_TEST_GENE.fasta")), False)

if __name__ == "__main__":
    unittest.main()