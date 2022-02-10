import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_extract_marker_genes():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/extract_marker_genes/data")
        expected_path = PurePosixPath(".tests/unit/extract_marker_genes/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        # Copy config file to temporary workdir
        sp.run([
            "cp",
            ".tests/config.yml",
            workdir
        ])
        # Copy test data to temporary workdir
        sp.run([
            "cp",
            ".tests/data/TEST",
            workdir
        ])
        # Copy run_assembly.txt to temporary workdir
        sp.run([
            "cp",
            ".tests/run_assembly.txt",
            workdir
        ])

        # dbg
        print("sequences/secE__GCF_000010525.1_ASM1052v1.fasta", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "sequences/secE__GCF_000010525.1_ASM1052v1.fasta",
            "-F", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        # Clean config, logs, run_assembly, and data from workdir
        sp.run([
            "rm",
            "-r",
            str(workdir) + "/config.yml",
            str(workdir) + "/logs",
            str(workdir) + "/TEST",
            str(workdir) + "/run_assembly.txt"
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
