#!/usr/bin/env python
import argparse
import json
from multiprocessing import Pool

import uproot
import numpy as np
import matplotlib.pyplot as plt
import opengate as gate

import itk
from itk import PCT as pct


def pv(verbose, *args, **kwargs):
    if verbose:
        print(*args, **kwargs)


def tof_mc(
    wepl,
    output,
    detector_distance_mm,
    number_of_particles,
    initial_energy,
    seed,
    verbose,
):

    u = gate.g4_units
    nm, mm, cm, m, sec, MeV = u.nm, u.mm, u.cm, u.m, u.second, u.MeV

    pv(verbose, "Starting simulation with following parameters: " + str(locals()))

    # Simulation
    sim = gate.Simulation()

    sim.random_engine = "MersenneTwister"
    sim.random_seed = "auto"
    sim.run_timing_intervals = [[0 * sec, 1 * sec]]
    sim.check_volumes_overlap = False
    sim.visu = False
    sim.visu_type = "vrml"
    sim.g4_verbose = False
    sim.progress_bar = verbose
    sim.number_of_threads = 1
    sim.random_seed = np.random.randint(65536) if seed is None else seed

    # Misc
    yellow = [1, 1, 0, 1]
    blue = [0, 0, 1, 1]

    # Geometry
    sim.world.material = "G4_AIR"
    sim.world.size = [4 * m, 4 * m, 4 * m]
    sim.world.color = [0, 0, 0, 0]

    # Phantom
    phantom_width_cm = 10.0
    if wepl > 0.0:
        phantom = sim.add_volume("Box", name="Phantom")
        phantom.size = [
            phantom_width_cm * cm,
            phantom_width_cm * cm,
            wepl * mm,
        ]
        phantom.material = "G4_WATER"
        phantom.color = blue
        phantom.set_max_step_size(1.0 * mm)

    # Beam
    source = sim.add_source("GenericSource", "mybeam")
    source.particle = "proton"
    source.energy.mono = initial_energy * MeV
    source.energy.type = "mono"
    source.position.type = "box"
    source.position.size = [1 * nm, 1 * nm, 1 * nm]
    source.position.translation = [
        0 * mm,
        0 * mm,
        (-detector_distance_mm / 2 - 10) * mm,
    ]
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, 1]
    source.n = number_of_particles

    # Physics list
    sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"
    sim.physics_manager.set_user_limits_particles(["proton"])

    # Phase spaces
    def add_detector(name, translation, attach_to_phantom=False):
        plane = sim.add_volume("Box", "PlanePhaseSpace" + name)
        plane.size = [2.0 * phantom_width_cm * cm, 2.0 * phantom_width_cm * cm, 1 * nm]
        plane.translation = translation
        plane.material = "G4_AIR"
        plane.color = yellow
        if attach_to_phantom:
            plane.mother = phantom.name

        phase_space = sim.add_actor("PhaseSpaceActor", "PhaseSpace" + name)
        phase_space.attached_to = plane.name
        phase_space.output_filename = (
            f"{output}/output/l{float(wepl):.3f}_ps{name}.root"
        )
        phase_space.attributes = [
            "LocalTime",
        ]
        if int(gate.utility.version("opengate").split(".")[1]) > 0:
            F = gate.actors.filters.GateFilterBuilder()
            phase_space.filter = F.ParticleName == "proton"
        else:
            particle_filter = sim.add_filter("ParticleFilter", "Filter" + name)
            particle_filter.particle = "proton"

    epsilon_mm = 1e-5
    add_detector("In", [0 * mm, 0 * mm, (-detector_distance_mm / 2 - epsilon_mm) * mm])
    add_detector("Out", [0 * mm, 0 * mm, (detector_distance_mm / 2 + epsilon_mm) * mm])

    sim.run()

    data_in = uproot.concatenate(
        f"{output}/output/l{float(wepl):.3f}_psIn.root", library="np"
    )
    data_out = uproot.concatenate(
        f"{output}/output/l{float(wepl):.3f}_psOut.root", library="np"
    )

    tofs = np.mean(data_out["LocalTime"]) - np.mean(data_in["LocalTime"])

    return tofs


def pctsdpweplfit(
    output,
    number_of_particles,
    wepl_samples,
    max_wepl,
    detector_distance,
    initial_energy,
    seed,
    verbose,
):
    wepls = np.linspace(0, max_wepl, wepl_samples, endpoint=False)

    results = []
    with Pool() as pool:
        for wepl in wepls:
            result = pool.apply_async(
                tof_mc,
                (
                    wepl,
                    output,
                    detector_distance,
                    number_of_particles,
                    initial_energy,
                    seed,
                    verbose,
                ),
            )
            results.append(result)
        pool.close()
        pool.join()
    tofs = [result.get() for result in results]

    dtofs = (tofs - tofs[0]) * 1000.0  # ps

    deg = 5
    coeffs = np.polyfit(wepls, dtofs, deg=deg)

    # Reproduce Figure 4 from Ulrich-Pur et al. [JINST, 2023]
    plt.figure()
    plt.scatter(wepls, dtofs, marker="+")
    xs = np.linspace(0.0, max_wepl, 1000)
    plt.plot(xs, np.polyval(coeffs, xs))
    plt.rcParams.update({"mathtext.default": "regular"})
    plt.xlabel("WEPL [mm]")
    plt.ylabel("$TOF-TOF_{air}$ [ps]")
    plt.savefig(f"{output}/fit.pdf")

    with open(f"{output}/fit.json", "w", encoding="utf-8") as f:
        json.dump(coeffs.tolist(), f)


def build_parser():
    parser = pct.PCTArgumentParser(
        description="Generate slowing-down power to WEPL calibration curve using the method from Ulrich-Pur et al. [JINST, 2023]",
    )

    parser.add_argument(
        "-o", "--output", help="Path of outputs", default="pctsdpweplfit"
    )
    parser.add_argument(
        "-n",
        "--number-of-particles",
        help="Number of generated particles",
        default=10000,
        type=int,
    )
    parser.add_argument(
        "--wepl-samples", help="Number of WEPL samples", type=int, default=20
    )
    parser.add_argument(
        "--max-wepl", help="Maximum WEPL to consider, in mm", type=int, default=200
    )
    parser.add_argument(
        "-d",
        "--detector-distance",
        help="Distance between detectors, in mm",
        default=150.0 + 2.0 * 100.0,  # Diameter of CTP404 + 2 times clearance
        type=float,
    )
    parser.add_argument(
        "-e",
        "--initial-energy",
        help="Initial energy of the protons, in MeV",
        default=200.0,
        type=float,
    )
    parser.add_argument("--seed", help="Seed for random number generator", type=int)
    parser.add_argument(
        "--verbose",
        "-v",
        help="Verbose execution",
        action="store_true",
    )

    return parser


def process(args_info: argparse.Namespace):
    pctsdpweplfit(
        output=args_info.output,
        number_of_particles=args_info.number_of_particles,
        wepl_samples=args_info.wepl_samples,
        max_wepl=args_info.max_wepl,
        detector_distance=args_info.detector_distance,
        initial_energy=args_info.initial_energy,
        seed=args_info.seed,
        verbose=args_info.verbose,
    )


def main(argv=None):
    parser = build_parser()
    args_info = parser.parse_args(argv)
    process(args_info)


if __name__ == "__main__":
    main()
