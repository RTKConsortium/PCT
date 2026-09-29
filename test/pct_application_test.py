import pytest
import itk
import urllib.request
import numpy as np
import uproot
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
        f"-o {output} --path-type phantom_length -d 220 -e 200 -l 220 --seed 1234"
    )

    tof_to_wepl_fit = np.loadtxt(output / "tof_to_wepl_fit_deg3.txt")
    reference = np.array(
        [
            8507.18830081492,
            -38342.09213481464,
            58268.25904916112,
            -29633.875319259325,
        ]
    )
    assert np.allclose(tof_to_wepl_fit, reference)
    eloss_to_wepl_fit = np.loadtxt(output / "eloss_to_wepl_fit_deg3.txt")
    reference = np.array(
        [
            -5.4569381663242e-06,
            -0.003455118013641362,
            2.236667157315751,
            -0.1768424416075149,
        ]
    )
    assert np.allclose(eloss_to_wepl_fit, reference)


def test_doublelut_application(tmp_path):
    output = tmp_path / "doublelut"
    pct.pctdoublelut(f"-o {output} -n 1000 --wepl-samples 10 -e 200 --seed 1234")

    tof_coeffs = np.loadtxt(output / "tof_coeffs_9.txt")
    tof_reference = np.array(
        [
            -1.304835138529633472e-19,
            1.356468922636278510e-16,
            -5.942720434270841880e-14,
            1.427297715745322610e-11,
            -2.045269455668975767e-09,
            1.780062155625161403e-07,
            -9.118495781246014652e-06,
            2.525195871208434648e-04,
            3.058213313702242298e-03,
            2.589942534059123773e-03,
        ]
    )
    assert np.allclose(tof_coeffs, tof_reference)

    vel_coeffs = np.loadtxt(output / "vel_coeffs_9.txt")
    vel_reference = np.array(
        [
            4.582905610876915593e-18,
            -4.742816184714826704e-15,
            2.062923179381136339e-12,
            -4.910573707663995126e-10,
            6.964091908878167104e-08,
            -5.995324018840303797e-06,
            3.031542041296401437e-04,
            -8.476856480304249125e-03,
            -4.918237198608615968e-02,
            1.696467164517849255e02,
        ]
    )
    assert np.allclose(vel_reference, vel_coeffs)


baseline_pairs_doublelut_mhd = download_file_fixture(
    "6a8f008e92f283f838800623", "baseline_pairs_doublelut.mhd"
)
baseline_pairs_doublelut_raw = download_file_fixture(
    "6a8f009092f283f838800626", "baseline_pairs_doublelut.raw"
)


def test_pairprotons_doublelut_application(
    tmp_path,
    phasespacein_root,
    phasespaceout_root,
    baseline_pairs_doublelut_mhd,
    baseline_pairs_doublelut_raw,
):
    output = tmp_path / "pairs_doublelut.mhd"

    tof_coeffs = tmp_path / "tof_coeffs.txt"
    np.savetxt(
        tof_coeffs,
        [
            2.796830534907338879e-22,
            -1.788065254786836053e-19,
            4.591994205164617155e-17,
            -5.662289174905150623e-15,
            3.899965359778268460e-13,
            -4.224885502649042288e-12,
            4.745768765394230579e-09,
            2.432353352721748280e-06,
            5.889704508365986406e-03,
            -9.553979783485672557e-07,
        ],
    )
    vel_coeffs = tmp_path / "vel_coeffs.txt"
    np.savetxt(
        vel_coeffs,
        [
            -1.811548969775848679e-19,
            1.487181098408996225e-16,
            -5.245527890416944471e-14,
            1.021392695091731880e-11,
            -1.214073609154529679e-09,
            8.787192875077365930e-08,
            -4.584117446224705820e-06,
            -1.447480766127036390e-04,
            -1.423014493400832636e-01,
            1.697332540903398126e02,
        ],
    )

    pct.pctpairprotons(
        f"-i {phasespacein_root} -j {phasespaceout_root} -o {output} --plane-in -110 --plane-out 110 --psin PhaseSpaceIn --psout PhaseSpaceOut --lut-tof {tof_coeffs} --lut-vel {vel_coeffs} --quadric 1 0 1 0 0 0 0 0 0 -10000 --angle 0"
    )

    output0000 = tmp_path / "pairs_doublelut0000.mhd"
    pairs_test = itk.array_from_image(itk.imread(output0000))
    pairs_baseline = itk.array_from_image(itk.imread(baseline_pairs_doublelut_mhd))
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
        f"-i {lomalinda_data} -o {output} --plane-in -167.2 --plane-out 167.2 --ps recoENTRY"
    )
    output0000 = str(output).replace(".", "0000.")
    test_lomalinda = itk.array_from_image(itk.imread(output0000))
    reference_lomalinda = itk.array_from_image(itk.imread(baseline_lomalinda_mhd))
    assert np.array_equal(test_lomalinda, reference_lomalinda)
    return output0000


baseline_addnoise = download_file_fixture(
    "6a8561052688ba21262c390a", "baseline_addnoise.root"
)


def test_addnoise_application(tmp_path, phasespacein_root, baseline_addnoise):
    output = tmp_path / "noise_test.root"
    tree = "PhaseSpaceIn"
    pct.pctaddnoise(
        f"-i {phasespacein_root} -o {output} --tree {tree} --material-budget .01 --tracker-distance 10 --translation -5 --noise-position 10 --noise-energy 10 --seed 1234"
    )
    root_test = uproot.open(output)[tree].arrays(library="np")
    root_baseline = uproot.open(baseline_addnoise)[tree].arrays(library="np")
    assert np.array_equal(root_test, root_baseline)
