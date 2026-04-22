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


lomalinda_data = download_file_fixture(
    "69e21803ed08a1c077afd077", "projection_045.root"
)
baseline_lomalinda_mhd = download_file_fixture(
    "69e89639ed08a1c077afd0d9", "baseline_lomalinda0000.mhd"
)
baseline_lomalinda_raw = download_file_fixture(
    "69e8963eed08a1c077afd0dc", "baseline_lomalinda0000.raw"
)


@pytest.fixture(scope="session")
def test_lomalinda_application(
    tmp_path_factory, lomalinda_data, baseline_lomalinda_mhd, baseline_lomalinda_raw
):
    output = tmp_path_factory.getbasetemp() / "lomalinda.mhd"
    pct.pctlomalinda(
        f"-i {lomalinda_data} -o {output} --plane-in -167.2 --plane-out 167.2 --ps recoENTRY -v"
    )
    output0000 = str(output).replace(".", "0000.")
    test_lomalinda = itk.array_from_image(itk.imread(output0000))
    reference_lomalinda = itk.array_from_image(itk.imread(baseline_lomalinda_mhd))
    assert np.array_equal(test_lomalinda, reference_lomalinda)
    return output0000


baseline_filterprotons_mhd = download_file_fixture(
    "69e89be4ed08a1c077afd0e5", "baseline_filterprotons.mhd"
)
baseline_filterprotons_raw = download_file_fixture(
    "69e89be5ed08a1c077afd0e8", "baseline_filterprotons.raw"
)


def test_filterprotons_application(
    tmp_path,
    test_lomalinda_application,
    baseline_filterprotons_mhd,
    baseline_filterprotons_raw,
):
    output = tmp_path / "filterprotons.mhd"
    pct.pctfilterprotons(
        f"-i {test_lomalinda_application} -o {output} --roi-radius 200 --fluence .5 -v --min-wepl 70 --max-wepl 180 --seed 1234"
    )
    test_filterprotons = itk.array_from_image(itk.imread(output))
    reference_filterprotons = itk.array_from_image(
        itk.imread(baseline_filterprotons_mhd)
    )
    assert np.array_equal(test_filterprotons, reference_filterprotons)
