import os
import shutil
import subprocess
import tempfile
import unittest

from genetools.command import main

class FilterGenesTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.index_fp = os.path.join(
            self.temp_dir, "TEST_TXIDS")
        self.index_contents = (
            "7\n"
            "9\n"
            "11\n"
            "14\n"
            "17\n"
            "19\n"
            "21\n"
            "23\n"
            "24\n"
            "25\n"
            "33\n")
        with open(self.index_fp, "w") as f:
            f.write(self.index_contents)
        
        self.gene_fp = os.path.join(
            self.temp_dir, "TEST_GENES")
        self.gene_contents = (
            "tRNA-Ser\n"
            "tRNA-Glu\n"
            "secE\n"
            "tsaE\n")
        with open(self.gene_fp, "w") as f:
            f.write(self.gene_contents)
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)

if __name__ == "__main__":
    unittest.main()