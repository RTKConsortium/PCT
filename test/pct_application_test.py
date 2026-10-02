import pytest
import itk
import urllib.request
import numpy as np
import uproot
from itk import PCT as pct


def download_file_fixture(file_hash, filename):

    @pytest.fixture(scope="session")
    def fixture(tmp_path_factory):
        path = tmp_path_factory.getbasetemp() / filename
        url = (
            f"https://data.kitware.com/api/v1/file/hashsum/sha512/{file_hash}/download"
        )
        with urllib.request.urlopen(url) as response, open(path, "wb") as out_file:
            out_file.write(response.read())
        return path

    return fixture


phasespacein_root = download_file_fixture(
    "59c88136fd0f6b00241fe0a32cde402f1843da4cbb7547a9ffdd1355f25155bb8a458dd041733448897440145e38c96a4facd0d976428a9e34561824cc34b2c6",
    "PhaseSpaceIn.root",
)
phasespaceout_root = download_file_fixture(
    "19d498c6d01bffac13b5aefe2d9382474ebdad77f64eeb536e6126a2ed2c296b02794d9ad78436d45745c51eeb85142cf87340e6967664c6851853063fde3ccd",
    "PhaseSpaceOut.root",
)
baseline_pairs_mhd = download_file_fixture(
    "5d3fdb78684c0355134bd2eb572f4303cab9e2848909cb4cd15e4b0fb42f461da9961c0523afd99c34d83ae211446aaa1baf11cad48c35f00abbeb481760ac20",
    "pairs0000.mhd",
)
baseline_pairs_raw = download_file_fixture(
    "e9eafd0490c52130452485a3b5328ac853f67996f13e58c444f4bfe745dd3b6f8ec371efbcffeb0b2570e303f9602c802e9fb496dce9388141762874c6029ab9",
    "pairs0000.raw",
)


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
    "90944573da740e633ffd1111994b4abd9f0542ae36fd2d760010d8942216bca0a3ec95c17f7500f7e9d30e26e5565d0ac8898d1de7ef5418366d5c82eb5b78c4",
    "baseline_pairs_doublelut.mhd",
)
baseline_pairs_doublelut_raw = download_file_fixture(
    "07322419bfc0998ae9744222a91c6f0d5b256b5bc1a1d8e6d6f9e5f7abd7682b56b9e22d9601cb261da41205601d24beb5d1e426f66a59fcd51f376637635103",
    "baseline_pairs_doublelut.raw",
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
    "fd8242eeddf2047f8e66ab357fb081df47e64c4e756358d3f15a9c7c7dfee32afa23a5873b5154226fa976e45c04b707356038e8474e9e4fc921332afdc43cd5",
    "projection_045.root",
)
baseline_lomalinda_mhd = download_file_fixture(
    "4a5439516da57ee7de5c482b2447ab43fa4eef318d97c667c6b6c8b959adaf9fc2a13ede4eb0ca3306eeedbf032bfb6f13cde79bb91b7e28676af987fe207f7f",
    "baseline_lomalinda0000.mhd",
)
baseline_lomalinda_raw = download_file_fixture(
    "bffdb9ff142eec5586509d575dc2c66b83689473bfb36a5c57299f91a5a46c2b032588e62647f528f550c4a13ec311f569582677a084f5c3f333412f25997a41",
    "baseline_lomalinda0000.raw",
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
    "20f798d1ddc57bd3f8784791cae050f665d47cd6c7aa097aee708a8ef007573e528cf9f84d90451dbdf3fe219d9adede2ae21cb36543a6924b786ab0c7325d31",
    "baseline_addnoise.root",
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
