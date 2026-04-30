import json
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


def test_weplfit_application(tmp_path):
    output = tmp_path / "weplfit"
    pct.pctweplfit(
        f"-o {output} --path-type phantom_length -d 220 -e 200 -l 220 --seed 1234  -v"
    )

    with open(output / "tof_to_wepl_fit_deg3.json", encoding="utf-8") as f:
        tof_to_wepl_fit = np.array(json.load(f))
        reference = np.array(
            [
                8507.18830081492,
                -38342.09213481464,
                58268.25904916112,
                -29633.875319259325,
            ]
        )
        assert np.allclose(tof_to_wepl_fit, reference)
    with open(output / "eloss_to_wepl_fit_deg3.json", encoding="utf-8") as f:
        eloss_to_wepl_fit = np.array(json.load(f))
        reference = np.array(
            [
                -5.4569381663242e-06,
                -0.003455118013641362,
                2.236667157315751,
                -0.1768424416075149,
            ]
        )
        assert np.allclose(eloss_to_wepl_fit, reference)
