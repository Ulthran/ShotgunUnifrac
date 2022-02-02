import os
import subprocess as sp

# If this test fails, install wget and/or gzip on your machine or
# use the containerized snakemake workflow
def test_wget_gzip():
    sp.check_output(["wget", "--help"])
    sp.check_output(["gzip", "--help"])