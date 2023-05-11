import os
import shutil
import tempfile
import unittest

from enum import Enum

from CorGE.command import main


class CommandTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = os.path.join(self.temp_dir, "output")

        # collect_genomes input
        self.data_dir = os.path.join("/".join(__file__.split("/")[:-1]), "test-data")
        self.ncbi_species_fp = os.path.join(self.data_dir, "TEST_TXIDS")
        self.ncbi_species_f = open(self.ncbi_species_fp)
        self.ncbi_accessions_fp = os.path.join(self.data_dir, "TEST_ACCS")
        self.ncbi_accessions_f = open(self.ncbi_accessions_fp)
        self.local_db_fp = os.path.join(self.data_dir, "TEST_LOCAL/")

        # collect_genomes outputs
        self.genomes_fp = os.path.join(self.temp_dir, "genomes/")

        # extract_genes input
        self.collected_genomes_fp = os.path.join(self.data_dir, "collected-genomes")

        # extract_genes output
        self.filtered_seqs_fp = os.path.join(self.temp_dir, "filtered-sequences")
        self.filtered_nucl_seqs_fp = os.path.join(
            self.temp_dir, "filtered-nucl-sequences"
        )
        self.merged_seqs_fp = os.path.join(self.temp_dir, "merged-sequences")

    def tearDown(self):
        self.ncbi_species_f.close()
        self.ncbi_accessions_f.close()

        shutil.rmtree(self.temp_dir)

    def test_collect_genomes(self):
        main(
            [
                "collect_genomes",
                "--output_fp",
                self.temp_dir,
                "--ncbi_species",
                self.ncbi_species_fp,
                "--ncbi_accessions",
                self.ncbi_accessions_fp,
                "--local",
                self.local_db_fp,
                # --outgroup left as default "2173"
            ]
        )

        self.assertEqual(
            os.listdir(self.genomes_fp).sort(),
            [
                "GCF_000012885.1.fna",
                "GCF_000007725.1.fna",
                "GCF_000020965.1.fna",
                "GCF_001735525.1.fna",
                "GCF_007197645.1.fna",
                "GCF_001375595.1.fna",
                "GCF_000218545.1.fna",
                "GCF_000010525.1.fna",
                "GCF_000378225.1.fna",
                "GCF_900111765.1.fna",
                "GCF_023159115.1.fna",
                "GCF_000016525.1.fna",
                "GCF_000012885.1.faa",
                "GCF_000007725.1.faa",
                "GCF_000020965.1.faa",
                "GCF_001735525.1.faa",
                "GCF_007197645.1.faa",
                "GCF_001375595.1.faa",
                "GCF_000218545.1.faa",
                "GCF_000010525.1.faa",
                "GCF_000378225.1.faa",
                "GCF_900111765.1.faa",
                "GCF_023159115.1.faa",
                "GCF_000016525.1.faa",
            ].sort(),
        )

    def test_extract_prot_genes(self):
        main(
            [
                "extract_genes",
                "--genomes",
                self.collected_genomes_fp,
                "--output",
                self.temp_dir,
            ]
        )

        import logging

        logging.warning(os.listdir(self.temp_dir))

        self.assertGreater(len(os.listdir(self.filtered_seqs_fp)), 700)
        self.assertEqual(len(os.listdir(self.merged_seqs_fp)), 71)

        with open(
            os.path.join(self.filtered_seqs_fp, "ADK__GCF_000007725.1.faa")
        ) as f:
            self.assertEqual(
                next(f).strip(),
                ">WP_011091539.1 nucleoside monophosphate kinase [Buchnera aphidicola]",
            )
            self.assertEqual(
                next(f).strip(),
                "MGFKCGFVGLPNVGKSTLFNYLTKLNIPADNYPFCTIKSNVGIVPVLDNRLNKIAQVVCSNKIIPATIELVDIAGLVKGAYKGEGLGNQFLDHIRDTNVIMHIVRCFENRYVTHIYGSVDPVRDVQIINLELILSDIEVCKNRMCKLEINKLSHNKQVNKELLILKKCVYHLEKSKSLRSLNLTEEEIFVINYLRLITLKPVVYIFNISIDQSRNLYKREIFDIIKNEHNAKTVNVCLDLMQSSKNDVSAYDHLSLKYKQLFNKMLKNVIWAGFNALNLITFFTAGKKEVHAWTTTNNLFIFQSVKCIHTDLSKGFIRAQVISYDDFIKYKGEKRSKELGKIRIEGKRYVICDGDIIHVLYNV",
            )

    def test_extract_nucl_genes(self):
        main(
            [
                "extract_genes",
                "--genomes",
                self.collected_genomes_fp,
                "--output",
                self.temp_dir,
                "--file_type",
                "nucl",
            ]
        )

        self.assertGreater(len(os.listdir(self.filtered_nucl_seqs_fp)), 700)
        self.assertEqual(len(os.listdir(self.merged_seqs_fp)), 71)

        with open(
            os.path.join(
                self.filtered_nucl_seqs_fp, "ADK__GCF_000007725.1.fna"
            )
        ) as f:
            self.assertEqual(
                next(f).strip(),
                ">lcl|NC_004545.1_cds_WP_011091539.1_407 [locus_tag=BBP_RS02150] [db_xref=GeneID:56470963] [protein=nucleoside monophosphate kinase] [protein_id=WP_011091539.1] [location=503441..504088] [gbkey=CDS]",
            )
            self.assertEqual(
                next(f).strip(),
                "ATGGGTTTTAAATGTGGTTTTGTTGGTTTACCTAATGTAGGAAAGTCTACTCTTTTTAATTATTTAACTAAATTAAATATTCCTGCAGATAATTATCCGTTTTGTACTATTAAGTCTAATGTTGGCATAGTTCCTGTTTTGGATAATCGCCTTAATAAAATAGCTCAAGTAGTTTGTTCTAATAAAATTATTCCAGCAACTATAGAGTTAGTAGATATTGCTGGATTAGTAAAAGGCGCTTATAAAGGTGAAGGATTAGGTAATCAATTTTTAGATCATATTAGAGACACAAATGTAATTATGCATATAGTGCGTTGTTTCGAGAATAGATATGTTACCCATATCTATGGTTCAGTAGATCCAGTGCGGGATGTACAAATTATAAATCTTGAGTTAATACTATCAGATATAGAAGTATGTAAAAATAGAATGTGCAAGCTTGAAATAAACAAGTTATCTCATAATAAACAAGTTAACAAGGAGTTATTAATATTAAAAAAATGCGTGTACCATTTAGAAAAAAGTAAAAGTTTGCGATCATTAAATTTAACTGAAGAAGAAATTTTTGTAATTAATTATTTAAGATTAATTACATTAAAACCTGTAGTGTATATTTTTAATATAAGTATAGATCAATCTAGAAATTTATATAAACGGGAAATTTTTGATATAATTAAAAATGAACATAATGCTAAAACAGTAAATGTTTGTTTAGATTTAATGCAAAGTAGCAAAAATGACGTTAGTGCATACGATCATCTTTCTTTAAAATATAAACAATTATTTAATAAAATGTTAAAAAATGTAATTTGGGCAGGTTTTAATGCTTTAAATTTAATTACTTTTTTTACTGCTGGCAAAAAAGAAGTTCATGCATGGACTACAACTAATAATTTGTTCATTTTTCAGTCTGTAAAATGTATTCATACAGATCTTAGTAAGGGATTCATTCGAGCTCAAGTGATTTCTTATGATGATTTTATTAAATATAAAGGAGAAAAGCGATCTAAAGAGTTAGGGAAAATTAGAATAGAAGGAAAAAGATATGTTATTTGTGATGGAGATATAATTCATGTTTTGTATAATGTATAG",
            )


if __name__ == "__main__":
    unittest.main()
