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
            ".tests/integration/full_run_snakemake/data/config.yml",
            os.path.join(tmpdir, "config.yml"),
        )
        shutil.copytree(
            "Astral/", os.path.join(tmpdir, "Astral/")
        )  # Would be better to have this defined in config
        shutil.copytree(
            ".tests/integration/full_run_snakemake/data/merged-sequences/",
            os.path.join(tmpdir, "merged-sequences/"),
        )

        with open(os.path.join(tmpdir, "config.yml"), "a") as f:
            f.write(f"\nDATA: {tmpdir}")

        os.system("conda config --set channel_priority strict")

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "all",
                "-c",
                "--conda-frontend",
                "conda",
                "--use-conda",
                "--conda-prefix",
                ".snakemake/",
                "--configfile",
                os.path.join(tmpdir, "config.yml"),
                # "--directory",
                # tmpdir,
            ]
        )

        # Check output
        if not os.path.exists(os.path.join(str(tmpdir), "RAxML_bestTree.outgroup")):
            raise ValueError("Full run did not produce RAxML_bestTree.outgroup")
