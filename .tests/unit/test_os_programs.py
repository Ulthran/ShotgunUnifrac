import os
import subprocess as sp

# If this test fails, install wget and/or gzip on your machine or
# use the containerized snakemake workflow
def test_os_programs():
    sp.check_output(["wget", "--help"])
    sp.check_output(["gzip", "--help"])
    sp.check_output(["cat", "--help"])