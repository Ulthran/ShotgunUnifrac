import os
import shutil
import subprocess
import tempfile
import unittest

from collections import Counter
from genetools.filterGenes import (get_txid, filter_seq_file, filter_seq_genes, extract_genes)

class FilterGenesTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        os.mkdir(os.path.join(self.temp_dir, "output"))
        os.mkdir(os.path.join(self.temp_dir, "output/ncbi"))
        os.mkdir(os.path.join(self.temp_dir, "output/sequences"))
        os.mkdir(os.path.join(self.temp_dir, "data"))

        self.run_assembly_fp = os.path.join(self.temp_dir, "data/run_assembly.txt")
        self.run_assembly_contents = (
            "assembly_accession	bioproject	biosample	wgs_master	refseq_category	taxid	species_taxid	organism_name	infraspecific_name	isolate	version_status	assembly_level	release_type	genome_rep	seq_rel_date	asm_name	submitter	gbrs_paired_asm	paired_asm_comp	ftp_path	excluded_from_refseq	relation_to_type_material	asm_not_live_date\n"
            "TEST_1	PRJNA224116	SAMN05384437	MCBT00000000.1	representative genome	23	23	Shewanella colwelliana	strain=CSB03KR		latest	Scaffold	Major	Full	2016/09/19	ASM173552v1	Chonnam National University	GCA_001735525.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/735/525/GCF_001735525.1_ASM173552v1			na\n"
            "TEST_2	PRJNA224116	SAMN05384437	MCBT00000000.1	representative genome	23	24	Shewanella colwelliana	strain=CSB03KR		latest	Scaffold	Major	Full	2016/09/19	ASM173552v1	Chonnam National University	GCA_001735525.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/735/525/GCF_001735525.1_ASM173552v1			na\n")
        with open(self.run_assembly_fp, "w") as f:
            f.write(self.run_assembly_contents)

        self.gene_file_fp = os.path.join(self.temp_dir, "data/TEST_GENES")
        self.gene_file_contents = (
            "TEST_a\n"
            "TEST_b\n")
        with open(self.gene_file_fp, "w") as f:
            f.write(self.gene_file_contents)

        # It's important that these test files be named something like *_*_cds/rna_from_genomic.fasta
        # to not through off parsing of the genome ID
        self.test_genome_cds_1_fp = os.path.join(self.temp_dir, "output/ncbi/TEST_1_cds_from_genomic.fasta")
        self.test_genome_rna_1_fp = os.path.join(self.temp_dir, "output/ncbi/TEST_1_rna_from_genomic.fasta")
        self.test_genome_cds_1_content = (
            ">TEST_A [gene=TEST_a]\n"
            "AAA\n"
            ">TEST_A [gene=DUMMY]\n"
            "GGG\n")
        self.test_genome_rna_1_content = (
            ">TEST_A [product=this produces TEST_b with magic]\n"
            "TTT\n"
            ">TEST_A [dummy=DUMMY] [protein=DUMMY]\n"
            "CCC\n")
        with open(self.test_genome_cds_1_fp, "w") as cds, open(self.test_genome_rna_1_fp, "w") as rna:
            cds.write(self.test_genome_cds_1_content)
            rna.write(self.test_genome_rna_1_content)

        self.test_genome_cds_2_fp = os.path.join(self.temp_dir, "output/ncbi/TEST_2_cds_from_genomic.fasta")
        self.test_genome_rna_2_fp = os.path.join(self.temp_dir, "output/ncbi/TEST_2_rna_from_genomic.fasta")
        self.test_genome_cds_2_content = (
            ">TEST_B [product=DUMMY] [protein=TEST_a] [gene=TEST_a]\n"
            "AAA\n"
            ">TEST_B [gene=DUMMY] [protein=DUMMY]\n"
            "GGG\n")
        self.test_genome_rna_2_content = (
            ">TEST_B [product=this produces DUMMY with magic]\n"
            "TTT\n"
            ">TEST_B [dummy=DUMMY] [protein=DUMMY]\n"
            "CCC\n")
        with open(self.test_genome_cds_2_fp, "w") as cds, open(self.test_genome_rna_2_fp, "w") as rna:
            cds.write(self.test_genome_cds_2_content)
            rna.write(self.test_genome_rna_2_content)
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_get_txid(self):
        self.assertEqual(get_txid("TEST_1", self.temp_dir), "23")
        self.assertEqual(get_txid("TEST_2", self.temp_dir), "24")

    def test_filter_seq_file(self):
        self.assertEqual(
            filter_seq_file(os.path.join(self.temp_dir, "output/ncbi/TEST_1_cds_from_genomic.fasta"), "TEST_1", "TEST_a"),
            [">23", "AAA"])
        self.assertEqual(
            filter_seq_file(os.path.join(self.temp_dir, "output/ncbi/TEST_1_rna_from_genomic.fasta"), "TEST_1", "TEST_b"),
            [">23", "TTT"])
        self.assertEqual(
            filter_seq_file(os.path.join(self.temp_dir, "output/ncbi/TEST_1_cds_from_genomic.fasta"), "TEST_1", "TEST_b"),
            None)
        
    def test_filter_seq_genes(self):
        self.assertEqual(
            filter_seq_genes("TEST_1", "TEST_a", os.path.join(self.temp_dir, "output/ncbi")),
            ['>23', 'AAA'])
        self.assertEqual(
            filter_seq_genes("TEST_1", "TEST_b", os.path.join(self.temp_dir, "output/ncbi")),
            ['>23', 'TTT'])
        self.assertEqual(
            filter_seq_genes("TEST_2", "TEST_a", os.path.join(self.temp_dir, "output/ncbi")),
            ['>24', 'AAA'])
        self.assertEqual(
            filter_seq_genes("TEST_2", "TEST_b", os.path.join(self.temp_dir, "output/ncbi")),
            [])
    
    def test_extract_genes(self):
        self.assertEqual(
            extract_genes(self.gene_file_fp, ["TEST_1", "TEST_2"], self.temp_dir),
            Counter({"TEST_a": 2, "TEST_b": 1}))

if __name__ == "__main__":
    unittest.main()