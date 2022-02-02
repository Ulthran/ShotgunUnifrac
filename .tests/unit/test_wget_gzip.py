import os
import subprocess as sp

def test_wget_gzip():
    sp.check_output(["wget", "--help"])
    sp.check_output(["gzip", "--help"])