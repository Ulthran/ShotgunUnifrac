import os
import shutil
import subprocess
import tempfile
import unittest

from CorGE.command import main

class CommandTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        # collect_genomes input
        self.data_dir = os.path.join('/'.join(__file__.split('/')[:-1]), "test-data")
        self.ncbi_species_fp = os.path.join(self.data_dir, "TEST_TXIDS")
        self.ncbi_accessions_fp = os.path.join(self.data_dir, "TEST_ACCS")
        self.local_db_fp = os.path.join(self.data_dir, "TEST_LOCAL/")

        # collect_genomes outputs
        self.nucl_fp = os.path.join(self.temp_dir, "nucleotide/")
        self.prot_fp = os.path.join(self.temp_dir, "protein/")
        self.outgroup_fp = os.path.join(self.temp_dir, "outgroup/")

        # extract_genes input
        self.collected_geneomes_fp = os.path.join(self.data_dir, "collected-genomes")

        # extract_genes output
        self.filtered_seqs_fp = os.path.join(self.temp_dir, "filtered-sequences")
        self.merged_seqs_fp = os.path.join(self.temp_dir, "merged-sequences")
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_collect_genomes(self):
        main([
            "collect_genomes",
            self.temp_dir,
            "--ncbi_species", self.ncbi_species_fp,
            "--ncbi_accessions", self.ncbi_accessions_fp,
            "--local", self.local_db_fp,
            # --outgroup left as default "2173"
        ])

        self.assertEqual(os.listdir(self.outgroup_fp).sort(), ['GCF_000016525.1.faa', 'GCF_000016525.1.fna'].sort())
        self.assertEqual(os.listdir(self.nucl_fp).sort(), ['GCF_000012885.1.fna', 'GCF_000007725.1.fna', 'GCF_000020965.1.fna',\
            'GCF_001735525.1.fna', 'GCF_007197645.1.fna', 'GCF_001375595.1.fna', 'GCF_000218545.1.fna', 'GCF_000010525.1\n.fna',\
            'GCF_000378225.1.fna', 'GCF_900111765.1.fna', 'GCF_023159115.1.fna'].sort())
        self.assertEqual(os.listdir(self.prot_fp).sort(), ['GCF_000012885.1.faa', 'GCF_000007725.1.faa', 'GCF_000020965.1.faa',\
            'GCF_001735525.1.faa', 'GCF_007197645.1.faa', 'GCF_001375595.1.faa', 'GCF_000218545.1.faa', 'GCF_000010525.1\n.faa',\
            'GCF_000378225.1.faa', 'GCF_900111765.1.faa', 'GCF_023159115.1.faa'].sort())

    def test_extract_genes(self):
        main([
            "extract_genes",
            self.collected_geneomes_fp,
            "--output", self.temp_dir,
        ])

        self.assertEqual(len(os.listdir(self.filtered_seqs_fp)), 488)
        self.assertEqual(len(os.listdir(self.merged_seqs_fp)), 41)

        with open(os.path.join(self.filtered_seqs_fp, "COG0012__GCF_000007725.1.faa")) as f:
            self.assertEqual(next(f).strip(), "> GCF_000007725.1")
            self.assertEqual(next(f).strip(), "MGFKCGFVGLPNVGKSTLFNYLTKLNIPADNYPFCTIKSNVGIVPVLDNRLNKIAQVVCSNKIIPATIELVDIAGLVKGA")


if __name__ == "__main__":
    unittest.main()