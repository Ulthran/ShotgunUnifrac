import os
import shutil
import sys

from pathlib import Path
import subprocess as sp
from tempfile import TemporaryDirectory

sys.path.insert(0, os.path.dirname(__file__))


def test_full_run():

    with TemporaryDirectory() as tmpdir_path:
        tmpdir = Path(tmpdir_path) / "workdir"
        os.makedirs(tmpdir)
        tmpdir = str(tmpdir)

        shutil.copyfile(
            ".tests/integration/full_supermatrix_run/data/config.yml",
            os.path.join(tmpdir, "config.yml"),
        )
        shutil.copytree(
            ".tests/integration/full_supermatrix_run/data/merged-sequences/",
            os.path.join(tmpdir, "merged-sequences/"),
        )

        os.system("conda config --set channel_priority strict")

        # Run the test job.
        output = sp.run(
            [
                "snakemake",
                "all",
                "-c",
                "--use-singularity",
                "--directory",
                tmpdir,
            ],
            capture_output=True,
        )

        # Check output
        print(output)
