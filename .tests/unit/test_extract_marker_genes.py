import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))
sys.path.append(str(os.path.dirname(__file__)) + "/../../workflow/scripts")

import common


def test_extract_marker_genes():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        sys.path.insert(0, str(workdir))
        data_path = PurePosixPath(".tests/unit/extract_marker_genes/data")
        expected_path = PurePosixPath(".tests/unit/extract_marker_genes/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Copy run_assembly.txt to temporary workdir (for matching genome id to taxon id)
        sp.run([
            "cp",
            ".tests/run_assembly.txt",
            str(workdir)
        ])

        # Write gene file
        with open(str(workdir) + "geneFile", "w") as geneF:
            geneF.write("secE")

        # Create downloaded genes
        os.system("mkdir " + str(workdir) + "/ncbi/")
        sp.run([
            "cp",
            "-r",
            ".tests/unit/extract_marker_genes/data/ncbi",
            str(workdir) + "/ncbi/"
        ])

        # Run the test job
        import downloadGenes as dg
        os.system("mkdir " + str(workdir) + "/sequences/")
        gene_counter = dg.extract_genes(str(workdir) + "geneFile", ["GCF_000010525.1_ASM1052v1"], str(workdir))

        # Check output
        if os.path.isfile(str(workdir) + "/sequences/secE__GCF_000010525.1_ASM1052v1.fasta"):
            if os.path.getsize(str(workdir) + "/sequences/secE__GCF_000010525.1_ASM1052v1.fasta") == 0:
                raise ValueError("Gene file is empty, extraction failed")
        else:
            raise ValueError("Gene file doesn't exist")

