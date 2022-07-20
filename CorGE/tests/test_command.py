import os
import shutil
import tempfile
import unittest

from enum import Enum

from CorGE.collect import collect_genomes
from CorGE.extract import extract_genes

class CommandTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        # collect_genomes input
        self.data_dir = os.path.join('/'.join(__file__.split('/')[:-1]), "test-data")
        self.ncbi_species_fp = os.path.join(self.data_dir, "TEST_TXIDS")
        self.ncbi_species_f = open(self.ncbi_species_fp)
        self.ncbi_accessions_fp = os.path.join(self.data_dir, "TEST_ACCS")
        self.ncbi_accessions_f = open(self.ncbi_accessions_fp)
        self.local_db_fp = os.path.join(self.data_dir, "TEST_LOCAL/")

        # collect_genomes outputs
        self.nucl_fp = os.path.join(self.temp_dir, "nucleotide/")
        self.prot_fp = os.path.join(self.temp_dir, "protein/")
        self.outgroup_fp = os.path.join(self.temp_dir, "outgroup/")

        # extract_genes input
        self.collected_genomes_fp = os.path.join(self.data_dir, "collected-genomes")

        # extract_genes output
        self.filtered_seqs_fp = os.path.join(self.temp_dir, "filtered-sequences")
        self.filtered_nucl_seqs_fp = os.path.join(self.temp_dir, "filtered-nucl-sequences")
        self.merged_seqs_fp = os.path.join(self.temp_dir, "merged-sequences")
    
    def tearDown(self):
        self.ncbi_species_f.close()
        self.ncbi_accessions_f.close()

        shutil.rmtree(self.temp_dir)
    
    def test_collect_genomes(self):
        # Doesn't seem to work with pytest-cov
        #main([
        #    "collect_genomes",
        #    self.temp_dir,
        #    "--ncbi_species", self.ncbi_species_fp,
        #    "--ncbi_accessions", self.ncbi_accessions_fp,
        #    "--local", self.local_db_fp,
        #    # --outgroup left as default "2173"
        #])

        collect_genomes(self.temp_dir, self.ncbi_species_f, self.ncbi_accessions_f, self.local_db_fp, "2173")

        self.assertEqual(os.listdir(self.outgroup_fp).sort(), ['GCF_000016525.1.faa', 'GCF_000016525.1.fna'].sort())
        self.assertEqual(os.listdir(self.nucl_fp).sort(), ['GCF_000012885.1.fna', 'GCF_000007725.1.fna', 'GCF_000020965.1.fna',\
            'GCF_001735525.1.fna', 'GCF_007197645.1.fna', 'GCF_001375595.1.fna', 'GCF_000218545.1.fna', 'GCF_000010525.1\n.fna',\
            'GCF_000378225.1.fna', 'GCF_900111765.1.fna', 'GCF_023159115.1.fna'].sort())
        self.assertEqual(os.listdir(self.prot_fp).sort(), ['GCF_000012885.1.faa', 'GCF_000007725.1.faa', 'GCF_000020965.1.faa',\
            'GCF_001735525.1.faa', 'GCF_007197645.1.faa', 'GCF_001375595.1.faa', 'GCF_000218545.1.faa', 'GCF_000010525.1\n.faa',\
            'GCF_000378225.1.faa', 'GCF_900111765.1.faa', 'GCF_023159115.1.faa'].sort())

    def test_extract_prot_genes(self):
        extract_genes(self.collected_genomes_fp, self.temp_dir, "prot")

        self.assertEqual(len(os.listdir(self.filtered_seqs_fp)), 488)
        self.assertEqual(len(os.listdir(self.merged_seqs_fp)), 41)

        with open(os.path.join(self.filtered_seqs_fp, "COG0012__GCF_000007725.1.faa")) as f:
            self.assertEqual(next(f).strip(), ">WP_011091313.1 redox-regulated ATPase YchF [Buchnera aphidicola]")
            self.assertEqual(next(f).strip(), "MGFKCGFVGLPNVGKSTLFNYLTKLNIPADNYPFCTIKSNVGIVPVLDNRLNKIAQVVCSNKIIPATIELVDIAGLVKGAYKGEGLGNQFLDHIRDTNVIMHIVRCFENRYVTHIYGSVDPVRDVQIINLELILSDIEVCKNRMCKLEINKLSHNKQVNKELLILKKCVYHLEKSKSLRSLNLTEEEIFVINYLRLITLKPVVYIFNISIDQSRNLYKREIFDIIKNEHNAKTVNVCLDLMQSSKNDVSAYDHLSLKYKQLFNKMLKNVIWAGFNALNLITFFTAGKKEVHAWTTTNNLFIFQSVKCIHTDLSKGFIRAQVISYDDFIKYKGEKRSKELGKIRIEGKRYVICDGDIIHVLYNV")

    def test_extract_nucl_genes(self):
        extract_genes(self.collected_genomes_fp, self.temp_dir, "nucl")

        self.assertEqual(len(os.listdir(self.filtered_nucl_seqs_fp)), 488)
        self.assertEqual(len(os.listdir(self.merged_seqs_fp)), 41)

        with open(os.path.join(self.filtered_nucl_seqs_fp, "COG0012__GCF_000007725.1.fna")) as f:
            self.assertEqual(next(f).strip(), ">lcl|NC_004545.1_cds_WP_011091313.1_169 [gene=ychF] [locus_tag=BBP_RS00905] [db_xref=GeneID:56470722] [protein=redox-regulated ATPase YchF] [protein_id=WP_011091313.1] [location=203069..204160] [gbkey=CDS]")
            self.assertEqual(next(f).strip(), "ATGGGTTTTAAATGTGGTTTTGTTGGTTTACCTAATGTAGGAAAGTCTACTCTTTTTAATTATTTAACTAAATTAAATATTCCTGCAGATAATTATCCGTTTTGTACTATTAAGTCTAATGTTGGCATAGTTCCTGTTTTGGATAATCGCCTTAATAAAATAGCTCAAGTAGTTTGTTCTAATAAAATTATTCCAGCAACTATAGAGTTAGTAGATATTGCTGGATTAGTAAAAGGCGCTTATAAAGGTGAAGGATTAGGTAATCAATTTTTAGATCATATTAGAGACACAAATGTAATTATGCATATAGTGCGTTGTTTCGAGAATAGATATGTTACCCATATCTATGGTTCAGTAGATCCAGTGCGGGATGTACAAATTATAAATCTTGAGTTAATACTATCAGATATAGAAGTATGTAAAAATAGAATGTGCAAGCTTGAAATAAACAAGTTATCTCATAATAAACAAGTTAACAAGGAGTTATTAATATTAAAAAAATGCGTGTACCATTTAGAAAAAAGTAAAAGTTTGCGATCATTAAATTTAACTGAAGAAGAAATTTTTGTAATTAATTATTTAAGATTAATTACATTAAAACCTGTAGTGTATATTTTTAATATAAGTATAGATCAATCTAGAAATTTATATAAACGGGAAATTTTTGATATAATTAAAAATGAACATAATGCTAAAACAGTAAATGTTTGTTTAGATTTAATGCAAAGTAGCAAAAATGACGTTAGTGCATACGATCATCTTTCTTTAAAATATAAACAATTATTTAATAAAATGTTAAAAAATGTAATTTGGGCAGGTTTTAATGCTTTAAATTTAATTACTTTTTTTACTGCTGGCAAAAAAGAAGTTCATGCATGGACTACAACTAATAATTTGTTCATTTTTCAGTCTGTAAAATGTATTCATACAGATCTTAGTAAGGGATTCATTCGAGCTCAAGTGATTTCTTATGATGATTTTATTAAATATAAAGGAGAAAAGCGATCTAAAGAGTTAGGGAAAATTAGAATAGAAGGAAAAAGATATGTTATTTGTGATGGAGATATAATTCATGTTTTGTATAATGTATAG")


if __name__ == "__main__":
    unittest.main()
