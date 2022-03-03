import os
import shutil
import subprocess
import tempfile
import unittest

from src.genetools.command import main

class CommandTests(unittest.TestCase):
    def setup(self):
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
        

    
    def teardown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_merge_and_config(self):
        main([
            self.index_fp,
            "--filter_genes", self.gene_fp,
            "--merge_genes",
            "--write_config",
            "--output_dir", self.temp_dir,
        ])

        with open(os.path.join(self.temp_dir, "output/config.yml")) as f:
            self.assertEqual(next(f), "# Config file for tree building pipeline")
            self.assertEqual(next(f), "# Genes to build trees from")
            self.assertEqual(next(f), "GENES: [\"tRNA-Ser\", \"tRNA-Glu\", \"secE\", \"tsaE\", ]")
        
        with open(os.path.join(self.temp_dir, "data/run_assembly.txt")) as f:
            self.assertEqual(next(f), "assembly_accession	bioproject	biosample	wgs_master	refseq_category	taxid	species_taxid	organism_name	infraspecific_name	isolate	version_status	assembly_level	release_type	genome_rep	seq_rel_date	asm_name	submitter	gbrs_paired_asm	paired_asm_comp	ftp_path	excluded_from_refseq	relation_to_type_material	asm_not_live_date")
            self.assertEqual(next(f), "GCF_001735525.1	PRJNA224116	SAMN05384437	MCBT00000000.1	representative genome	23	23	Shewanella colwelliana	strain=CSB03KR		latest	Scaffold	Major	Full	2016/09/19	ASM173552v1	Chonnam National University	GCA_001735525.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/735/525/GCF_001735525.1_ASM173552v1			na")

        self.assertEqual(os.path.exists(os.path.join(self.temp_dir, "output/ncbi/GCF_001735525.1_ASM173552v1.fasta")), True)


if __name__ == "__main__":
    unittest.main()