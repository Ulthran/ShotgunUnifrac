import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))
sys.path.append(str(os.path.dirname(__file__)) + "/../../workflow/scripts")

from collections import Counter


def test_merge_marker_genes():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/merge_marker_genes/data")
        expected_path = PurePosixPath(".tests/unit/merge_marker_genes/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        
        # write gene extraction file to be "merged" into gene file
        os.system("mkdir " + str(workdir) + "/sequences/")
        os.system("mkdir " + str(workdir) + "/merged-sequences/")
        with open(str(workdir) + "/sequences/secE.fasta", "w") as secE:
            secE.write(">lcl|NC_009937.1_cds_WP_043878890.1_905 [gene=secE] [locus_tag=AZC_RS04525] [protein=preprotein translocase subunit SecE] [protein_id=WP_043878890.1] [location=983866..984063] [gbkey=CDS]")
            secE.write("ATGGCAAAGAACAGTCCCGTGGAGTTCTTCCAGCAGGTCCGCACCGAGACGGCGAAGGTGACCTGGCCGTCCCGGCGCGAGACGCTGATCACCACCGCCATGGTCTTCGTCATGGTGCTGCTGGCGTCCATCTTCTTCCTGGTCGTGGACCAGATCCTGCGATTCGGCGTCAGCCAGATCCTCAGCATCGGCCATTGA")

        # Run the test job
        import downloadGenes as dg
        gene_list = dg.merge_genes(Counter({"secE": 4}), str(workdir))

        if not gene_list:
            raise ValueError("gene_list empty, something went wrong in merge_genes function")
        if not os.path.isfile(str(workdir) + "/merged-sequences/secE.fasta"):
            raise ValueError("Failed to create merged gene file")
