import os
import shutil
import subprocess
import tempfile
import unittest

from genetools.downloadGenes import (check_for_assembly, download_assembly, find_genome_accessions, prepare_run_assembly, download_genomes)

class DownloadGenesTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        os.mkdir(os.path.join(self.temp_dir, "data"))
        os.mkdir(os.path.join(self.temp_dir, "output"))
        os.mkdir(os.path.join(self.temp_dir, "output/ncbi"))

        self.run_assembly_fp = os.path.join(self.temp_dir, "data/run_assembly.txt")
        self.log_f = open(os.path.join(self.temp_dir, "log"), "w")

        self.taxon_list_fp = os.path.join(self.temp_dir, "data/TEST_TXIDS")
        self.taxon_list_content = (
            "# TITLE LINE\n"
            "23\n"
            "24\n"
            "NOT_FOUND\n")
        with open(self.taxon_list_fp, "w") as f:
            f.write(self.taxon_list_content)
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_find_genome_accessions(self):
        None if check_for_assembly(self.temp_dir) else download_assembly(self.temp_dir)
        self.assertEqual(
            find_genome_accessions(["23", "24", "NOT_FOUND"], self.temp_dir),
            ["NOT_FOUND"])
    
    def test_prepare_run_assembly(self):
        failed_txids = prepare_run_assembly(self.taxon_list_fp, self.temp_dir, self.log_f)
        self.assertEqual(failed_txids, ["NOT_FOUND"])
        with open(self.run_assembly_fp) as f:
            self.assertEqual(
                "assembly_accession	bioproject	biosample	wgs_master	refseq_category	taxid	species_taxid	organism_name	infraspecific_name	isolate	version_status	assembly_level	release_type	genome_rep	seq_rel_date	asm_name	submitter	gbrs_paired_asm	paired_asm_comp	ftp_path	excluded_from_refseq	relation_to_type_material	asm_not_live_date\n",
                next(f))
            self.assertEqual(
                "GCF_001735525.1	PRJNA224116	SAMN05384437	MCBT00000000.1	representative genome	23	23	Shewanella colwelliana	strain=CSB03KR		latest	Scaffold	Major	Full	2016/09/19	ASM173552v1	Chonnam National University	GCA_001735525.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/735/525/GCF_001735525.1_ASM173552v1			na\n",
                next(f))
            self.assertEqual(
                "GCF_016406325.1	PRJNA224116	SAMN17120581		representative genome	24	24	Shewanella putrefaciens	strain=ATCC 8071		latest	Complete Genome	Major	Full	2020/12/27	ASM1640632v1	Ocean University of China	GCA_016406325.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/406/325/GCF_016406325.1_ASM1640632v1		assembly from type material	na\n",
                next(f))
    
    def test_download_genomes(self):
        with open(self.run_assembly_fp, "w") as f:
            f.write("assembly_accession	bioproject	biosample	wgs_master	refseq_category	taxid	species_taxid	organism_name	infraspecific_name	isolate	version_status	assembly_level	release_type	genome_rep	seq_rel_date	asm_name	submitter	gbrs_paired_asm	paired_asm_comp	ftp_path	excluded_from_refseq	relation_to_type_material	asm_not_live_date\n")
            f.write("GCF_001735525.1	PRJNA224116	SAMN05384437	MCBT00000000.1	representative genome	23	23	Shewanella colwelliana	strain=CSB03KR		latest	Scaffold	Major	Full	2016/09/19	ASM173552v1	Chonnam National University	GCA_001735525.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/735/525/GCF_001735525.1_ASM173552v1			na\n")
            f.write("GCF_016406325.1	PRJNA224116	SAMN17120581		representative genome	24	24	Shewanella putrefaciens	strain=ATCC 8071		latest	Complete Genome	Major	Full	2020/12/27	ASM1640632v1	Ocean University of China	GCA_016406325.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/406/325/GCF_016406325.1_ASM1640632v1		assembly from type material	na\n")
            f.write("FAILED_DOWNLOAD	PRJNA224116	SAMN17120581		representative genome	24	0	Shewanella putrefaciens	strain=ATCC 8071		latest	Complete Genome	Major	Full	2020/12/27	ASM1640632v1	Ocean University of China	GCA_016406325.1	identical	https://badurl.stinky		assembly from type material	na\n")

        downloaded_genome_ids, failed_genome_ids = download_genomes(self.temp_dir, self.log_f)
        self.assertEqual(downloaded_genome_ids, ["GCF_001735525.1_ASM173552v1", "GCF_016406325.1_ASM1640632v1"])
        self.assertEqual(failed_genome_ids, ["badurl.stinky"])

        self.assertTrue(os.path.exists(os.path.join(self.temp_dir, "output/ncbi/GCF_001735525.1_ASM173552v1_cds_from_genomic.fasta")))
        self.assertTrue(os.path.exists(os.path.join(self.temp_dir, "output/ncbi/GCF_001735525.1_ASM173552v1_rna_from_genomic.fasta")))
        self.assertTrue(os.path.exists(os.path.join(self.temp_dir, "output/ncbi/GCF_016406325.1_ASM1640632v1_cds_from_genomic.fasta")))
        self.assertTrue(os.path.exists(os.path.join(self.temp_dir, "output/ncbi/GCF_016406325.1_ASM1640632v1_rna_from_genomic.fasta")))

if __name__ == "__main__":
    unittest.main()