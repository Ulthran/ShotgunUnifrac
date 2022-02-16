import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))


def test_full_run():

    with TemporaryDirectory() as tmpdir:
        # Create temporary workdir
        workdir = Path(tmpdir) / "workdir"
        sp.run([
            "mkdir",
            workdir
        ])
        # Copy config file to temporary workdir
        sp.run([
            "cp",
            ".tests/config.yml",
            workdir
        ])
        # Copy run_assembly.txt to temporary workdir
        sp.run([
            "cp",
            ".tests/run_assembly.txt",
            workdir
        ])
	    # Copy Astral to temporary workdir
        sp.run([
            "cp",
            "-r",
            "Astral",
            workdir
        ])

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "-c",
            "--use-conda",
            "--directory",
            workdir,
        ])

        # Clean config, logs, run_assembly, and data from workdir
        sp.run([
            "rm",
            "-r",
            str(workdir) + "/config.yml",
            str(workdir) + "/logs",
            str(workdir) + "/run_assembly.txt",
            str(workdir) + "/Astral"
        ])

        # Check output
        if not os.path.exists(str(workdir) + "/RAxML_rootedTree.final"):
            raise ValueError("Full run did not produce RAxML_rootedTree.final")
