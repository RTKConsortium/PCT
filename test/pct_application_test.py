import pytest
import itk
import urllib.request
import numpy as np
from itk import PCT as pct


def download_file_fixture(file_key, filename):

    @pytest.fixture(scope="session")
    def fixture(tmp_path_factory):
        path = tmp_path_factory.getbasetemp() / filename
        url = f"https://data.kitware.com/api/v1/file/{file_key}/download"
        with urllib.request.urlopen(url) as response, open(path, "wb") as out_file:
            out_file.write(response.read())
        return path

    return fixture


phasespacein_root = download_file_fixture(
    "69cbef28303cec2e64feb78a", "PhaseSpaceIn.root"
)
phasespaceout_root = download_file_fixture(
    "69cbef2a303cec2e64feb78d", "PhaseSpaceOut.root"
)
baseline_pairs_mhd = download_file_fixture("69cbefdf303cec2e64feb790", "pairs0000.mhd")
baseline_pairs_raw = download_file_fixture("69cbefe0303cec2e64feb793", "pairs0000.raw")


def test_pairprotons_application(
    tmp_path,
    phasespacein_root,
    phasespaceout_root,
    baseline_pairs_mhd,
    baseline_pairs_raw,
):
    output = tmp_path / "pairs_test.mhd"
    pct.pctpairprotons(
        f"-i {phasespacein_root} -j {phasespaceout_root} -o {output} --plane-in -110 --plane-out 110 --psin PhaseSpaceIn --psout PhaseSpaceOut"
    )
    output0000 = tmp_path / "pairs_test0000.mhd"
    pairs_test = itk.array_from_image(itk.imread(output0000))
    pairs_baseline = itk.array_from_image(itk.imread(baseline_pairs_mhd))
    assert np.array_equal(pairs_test, pairs_baseline)
