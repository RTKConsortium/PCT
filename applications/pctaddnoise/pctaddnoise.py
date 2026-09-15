"""
Add uncertainties to proton CT tracker data.

Proton CT detectors are typically simulated in GATE using a single PhaseSpaceActor. However, proton CT trackers actually consist of two detector planes. The position of the proton is detected twice, which allows to additionally detect the proton direction. However, several effects affect the measured positions and directions: protons may undergo scattering in the detector, and proton CT trackers are usually made out of strips that discretize the detection area.

This application takes as input a ROOT file generated from a PhaseSpaceActor representing a proton CT tracker. The following uncertainties are taken into account:
- if the ROOT file contains position and direction branches: realistic position and direction uncertainties based on Krah et al. (PMB, 2018)
- if the ROOT file contains an energy branch: Gaussian energy uncertainty
- if the ROOT file contains a time branch: Gaussian time uncertainty

This script is adapted from "AddTrackerUncertainty.py" (https://github.com/RTKConsortium/PCT/blob/6afe8ee0a0c25fd3e892761a237f120fd19198e4/AddTrackerUncertainty.py) from earlier versions of PCT.
"""

#!/usr/bin/env python
import sys
import numpy as np
import hepunits
import uproot
import argparse
import itk
from itk import PCT as pct


def build_parser():
    parser = pct.PCTArgumentParser(
        description="Add uncertainties to proton CT tracker data"
    )
    parser.add_argument(
        "--material-budget", default=5e-3, type=float, help="Material budget (unitless)"
    )
    parser.add_argument(
        "--tracker-distance",
        default=10.0,
        type=float,
        help="Distance between trackers (cm)",
    )
    parser.add_argument("-i", "--input", required=True, help="Input file name")
    parser.add_argument("--tree", required=True, help="Name of tree in ROOT file")
    parser.add_argument("-o", "--output", required=True, help="Output file name")
    parser.add_argument(
        "--translation",
        default=0.0,
        help="Translation of the detector position in the direction perpendicular to the detector (mm)",
    )
    parser.add_argument(
        "--noise-position",
        help="Standard deviation of the Gaussian noise on the position",
        type=float,
    )
    parser.add_argument(
        "--noise-energy",
        help="Standard deviation of the Gaussian noise on the energy",
        type=float,
    )
    parser.add_argument(
        "--noise-time",
        help="Standard deviation of the Gaussian noise on the time",
        type=float,
    )
    parser.add_argument("--seed", help="Random seed", type=int)
    parser.add_argument(
        "--verbose", "-v", help="Verbose execution", default=False, action="store_true"
    )
    return parser


def get_sigma_sc(energy, x_over_x0, sp, dt):
    """
    Get the Σ_sc matrix as defined by Equation (27) of Krah et al. (PMB, 2018).

    Args:
        - energy: energies of the protons (MeV).
        - x_over_x0: material budget.
        - sp: standard deviation of the tracker uncertainty (mm).
        - dt: distance between trackers (cm).
    """
    if sp is None:
        return

    proton_mass_c2 = 938.272013 * hepunits.MeV
    betap = (energy + 2 * proton_mass_c2) * energy / (energy + proton_mass_c2)

    # Equation (25) and (26) from Krah et al. (PMB, 2018).
    T = np.zeros((2, 2))
    T[0, 1] = 1
    T[1, 0] = -1 / dt
    T[1, 1] = 1 / dt

    sigma = (
        13.6
        * hepunits.MeV
        / betap
        * np.sqrt(x_over_x0)
        * (1 + 0.038 * np.log(x_over_x0))
    )
    sigma_sc = np.zeros((energy.size, 2, 2))
    sigma_sc[:, 1, 1] = sigma**2
    return np.tile(sp**2 * T @ T.T, (energy.size, 1, 1)) + sigma_sc


def add_tracker_uncertainty(
    data, rng, material_budget, noise_position, tracker_distance
):
    """
    Add tracker uncertainties as detailed in Krah et al. (PMB, 2018), section 2.5.

    Args:
        - data: data from the ROOT file, in NumPy format.
        - rng: random number generator.
        - material_budget: material budget.
        - noise_position: standard deviation of the tracker uncertainty (mm).
        - tracker_distance: distance between trackers (cm).
    """
    e = data["KineticEnergy"]
    sigma = get_sigma_sc(
        e, material_budget, noise_position * hepunits.mm, tracker_distance * hepunits.cm
    )
    w, q = np.linalg.eig(np.linalg.inv(sigma))
    q = np.real(q)
    xr = rng.standard_normal((e.size, 2, 2))
    W = np.zeros((e.size, 2, 2))
    W[:, 0, 0] = 1.0 / np.sqrt(w[:, 0])
    W[:, 1, 1] = 1.0 / np.sqrt(w[:, 1])
    dy_uncert = np.matmul(np.matmul(q, W), xr)
    data["Position_X"] += dy_uncert[:, 0, 0]
    data["Position_Y"] += dy_uncert[:, 0, 1]
    data["Direction_X"] += dy_uncert[:, 1, 0]
    data["Direction_Y"] += dy_uncert[:, 1, 1]


def add_gaussian_noise(data, branch, rng, noise, clamp=None):
    """
    Add Gaussian noise to an arbitrary branch of ROOT data.

    Args:
        - data: data from the ROOT file, in NumPy format.
        - branch: what branch to consider.
        - rng: random number generator.
        - noise: standard deviation.
        - clamp: values below this number will be clamped to this number.
    """
    if noise is None:
        return

    try:
        data[branch] += rng.normal(scale=noise, size=len(data["KineticEnergy"]))
        if clamp is not None:
            data[branch] = np.where(data[branch] < clamp, clamp, data[branch])
    except KeyError:
        print(
            f"Warning: cannot apply noise of {noise} on branch {branch} as the branch does not exist in the ROOT file! Skipping.",
            file=sys.stderr,
        )


def process(args_info: argparse.Namespace):

    rng = np.random.default_rng(args_info.seed)

    if args_info.verbose:
        print("Reading input ROOT file")
    data = uproot.open(args_info.input)[args_info.tree].arrays(library="np")

    # Move to entrance and exit detector (new) positions
    if args_info.verbose and args_info.translation is not None:
        print("Applying translations…")
        for pos in ["X", "Y", "Z"]:
            data[f"Position_{pos}"] += (
                args_info.translation / data["Direction_Z"]
            ) * data[f"Direction_{pos}"]

    if args_info.noise_position is not None and args_info.noise_position > 0.0:
        if args_info.verbose:
            print("Applying noise…")
        add_tracker_uncertainty(
            data,
            rng,
            args_info.material_budget,
            args_info.noise_position,
            args_info.tracker_distance,
        )

    add_gaussian_noise(data, "KineticEnergy", rng, args_info.noise_energy, clamp=0.0)
    add_gaussian_noise(data, "GlobalTime", rng, args_info.noise_time)

    if args_info.verbose:
        print("Writing output ROOT file…")
    with uproot.recreate(args_info.output) as output_file:
        output_file[args_info.tree] = data


def main(argv=None):
    parser = build_parser()
    args_info = parser.parse_args(argv)
    process(args_info)


if __name__ == "__main__":
    main()
