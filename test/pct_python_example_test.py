import sys
import os
import pytest
import urllib.request
import runpy
from pathlib import Path

# Base examples directory
EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def run_example(tmp_path, rel_script, *args):
    script = EXAMPLES / rel_script
    os.chdir(tmp_path)
    sys.argv = [str(script), *map(str, args)]
    runpy.run_path(str(script), run_name="__main__")


def test_Reconstruction(tmp_path):
    run_example(tmp_path, "Reconstruction/Reconstruction.py", "Reconstruction")
