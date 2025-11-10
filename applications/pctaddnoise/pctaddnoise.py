#!/usr/bin/env python
import argparse
import json
import sys
import itk
from itk import PCT as pct
import numpy as np
import numpy.lib.recfunctions as rfn
import hepunits


def build_parser():
    parser = pct.PCTArgumentParser(description="Add noise to ROOT data")
    parser.add_argument(
        "-i",
        "--input",
        help="Root phase space file of particles",
        required=True,
    )
    parser.add_argument("-t", "--tree", help="Root tree name", required=True)
    parser.add_argument("-o", "--output", help="Output file name", required=True)

    parser.add_argument("--seed", help="Random number generator seed", type=int)
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

    parser.add_argument(
        "--verbose", "-v", help="Verbose execution", default=False, action="store_true"
    )

    tracker_uncertainty_group = parser.add_argument_group("Tracker uncertainty")
    tracker_uncertainty_group.add_argument(
        "--material-budget", help="Material budget x/x0", type=float, default=5e-3
    )
    tracker_uncertainty_group.add_argument(
        "--noise-position",
        help="Standard deviation of the tracker uncertainty (mm)",
        type=float,
    )
    tracker_uncertainty_group.add_argument(
        "--tracker-distance", help="Distance between trackers (cm)", type=float
    )

    return parser


def get_sigma_sc(energy, x_over_x0, sp, dt):
    proton_mass_c2 = 938.272013 * hepunits.MeV
    betap = (energy + 2 * proton_mass_c2) * energy / (energy + proton_mass_c2)

    # Equation (25) and (26) from [Krah et al, PMB, 2018]:
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
    e = data["KineticEnergy"]
    sigma = get_sigma_sc(
        e, material_budget, noise_position * hepunits.mm, tracker_distance * hepunits.cm
    )
    w, q = np.linalg.eig(np.linalg.inv(sigma))
    xr = rng.standard_normal((e.size, 2, 2))
    W = np.zeros((e.size, 2, 2))
    W[:, 0, 0] = 1.0 / np.sqrt(w[:, 0])
    W[:, 1, 1] = 1.0 / np.sqrt(w[:, 1])
    dy_uncert = np.matmul(np.matmul(q, W), xr)
    data["Position_X"] += dy_uncert[:, 0, 0]
    data["Position_Y"] += dy_uncert[:, 0, 1]
    data["Direction_X"] += dy_uncert[:, 1, 0]
    data["Direction_Y"] += dy_uncert[:, 1, 1]


def add_gaussian_noise(data, branch, rng, noise):
    if noise is not None:
        try:
            data[branch] += rng.normal(scale=noise, size=len(data["KineticEnergy"]))
        except KeyError:
            print(
                f"Warning: cannot apply noise of {noise} on branch {branch} as the branch does not exist in the Root file! Skipping.",
                file=sys.stderr,
            )


def process(args_info: argparse.Namespace):
    import uproot

    rng = np.random.default_rng(args_info.seed)

    if args_info.verbose:
        print("Reading input Root file")
    tree = uproot.open(args_info.input)[args_info.tree]
    data = tree.arrays(library="np")

    if args_info.verbose:
        print("Applying noise…")

    add_tracker_uncertainty(
        data,
        rng,
        args_info.material_budget,
        args_info.noise_position,
        args_info.tracker_distance,
    )

    add_gaussian_noise(data, "KineticEnergy", rng, args_info.noise_energy)
    add_gaussian_noise(data, "LocalTime", rng, args_info.noise_time)

    if args_info.verbose:
        print("Writing output Root file…")
    with uproot.recreate(args_info.output) as output_file:
        output_file[args_info.tree] = data


def main(argv=None):
    parser = build_parser()
    args_info = parser.parse_args(argv)
    process(args_info)


if __name__ == "__main__":
    main()
