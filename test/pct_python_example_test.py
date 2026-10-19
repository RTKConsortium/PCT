import sys
import os
import pytest
import urllib.request
import runpy
from pathlib import Path
import warnings

# Base examples directory
EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def run_example(tmp_path, rel_script, *args):
    script = EXAMPLES / rel_script
    os.chdir(tmp_path)
    sys.argv = [str(script), *map(str, args)]
    runpy.run_path(str(script), run_name="__main__")


def test_GradientDescentReconstruction(tmp_path):
    # Catch DeprecationWarnings coming from opengate and making PCT tests fail
    # To be removed once the warnings are gone from opengate
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        run_example(
            tmp_path,
            "GradientDescentReconstruction/GradientDescentReconstruction.py",
            "GradientDescentReconstruction",
        )
