#!/usr/bin/env python
import argparse
import sys
from multiprocessing import Pool, Manager
import os
import warnings

import numpy as np
import itk
from itk import PCT as pct


def build_parser():
    parser = pct.PCTArgumentParser(
        description="Convert TOF to WEPL using a fit on Monte Carlo data",
    )

    parser.add_argument(
        "-o", "--output", help="Path of outputs", default="pctdoublelut"
    )
    parser.add_argument(
        "-n",
        "--number-of-particles",
        help="Number of generated particles",
        default=1000,
        type=int,
    )
    parser.add_argument("-m", "--material", help="Material", default="G4_WATER")
    parser.add_argument(
        "--wepl-samples",
        help="Number of WEPL samples",
        default=(2 * 260) + 1,
        type=int,
    )
    parser.add_argument(
        "--max-wepl",
        help="Maximum WEPL",
        default=260.0,
        type=float,
    )
    parser.add_argument(
        "-w",
        "--phantom-width",
        help="Phantom width",
        default=40.0,
        type=float,
    )
    parser.add_argument(
        "-e",
        "--initial-energy",
        help="Initial energy of the protons (in MeV)",
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


epsilon_mm = 1e-5


def pv(verbose, *args, **kwargs):
    if verbose:
        print(*args, **kwargs)


def get_gate_output(output):
    return os.path.join(output, "gate")


def tof_fit_mc(
    wepl,
    output,
    material,
    phantom_width_cm,
    number_of_particles,
    initial_energy,
    seed,
    verbose,
):
    import opengate as gate

    u = gate.g4_units
    nm, mm, cm, m, sec, MeV = u.nm, u.mm, u.cm, u.m, u.second, u.MeV

    # Simulation
    sim = gate.Simulation()

    sim.random_engine = "MersenneTwister"
    sim.random_seed = seed
    sim.run_timing_intervals = [[0.0 * sec, 1.0 * sec]]
    sim.check_volumes_overlap = False
    sim.g4_verbose = True
    sim.progress_bar = verbose

    # Misc
    yellow = [1, 1, 0, 1]
    blue = [0, 0, 1, 1]

    # Geometry
    sim.world.material = "G4_AIR"
    sim.world.size = [4 * m, 4 * m, 4 * m]

    # Phantom
    if wepl > 0.0:
        phantom = sim.add_volume("Box", name="Phantom")
        phantom.size = [
            phantom_width_cm * cm,
            phantom_width_cm * cm,
            wepl * mm,
        ]
        phantom.material = material
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
        (-wepl / 2 - 0.1) * mm,
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
        plane.size = [phantom_width_cm * cm, phantom_width_cm * cm, 1 * nm]
        plane.translation = translation
        plane.material = "G4_AIR"
        plane.color = yellow
        if attach_to_phantom:
            plane.mother = phantom.name

        phase_space = sim.add_actor("PhaseSpaceActor", "PhaseSpace" + name)
        phase_space.attached_to = plane.name
        phase_space.output_filename = os.path.join(
            get_gate_output(output), f"wepl{wepl:.3f}_ps{name}.root"
        )
        phase_space.attributes = [
            "EventID",
            "TrackID",
            "Position",
            "LocalTime",
            "PostVelocity",
        ]
        if int(gate.utility.version("opengate").split(".")[1]) > 0:
            F = gate.actors.filters.GateFilterBuilder()
            phase_space.filter = F.ParticleName == "proton"
        else:
            particle_filter = sim.add_filter("ParticleFilter", "Filter" + name)
            particle_filter.particle = "proton"

    add_detector("In", [0.0 * mm, 0.0 * mm, (-wepl / 2) * mm])
    add_detector("Out", [0.0 * mm, 0.0 * mm, (wepl / 2) * mm])

    sim.run()


def process_wepl(wepl, output, verbose):
    import uproot

    wepls = []
    tofs_wepl = []
    velocities_wepl = []

    data = uproot.concatenate(
        os.path.join(get_gate_output(output), f"wepl{wepl:.3f}_*.root"), library="np"
    )
    pv(
        verbose,
        "Loaded",
        len(data["EventID"]),
        "events for WEPL",
        wepl,
    )

    # Sort data if needed
    ws = data["Position_Z"]
    if not np.all(ws[1:] > ws[:-1]):
        pv(verbose, "Sorting input data…")
        index_sorted = np.argsort(ws)
        for key in data.keys():
            data[key] = data[key][index_sorted]

    for n in np.unique(data["EventID"]):
        event_mask = data["EventID"] == n

        number_of_hits = np.sum(event_mask)
        if number_of_hits == 0:
            continue

        times = data["LocalTime"][event_mask]
        tof = times[-1] - times[0]

        velocity = data["PostVelocity"][event_mask][-1]

        wepls.append(wepl)
        tofs_wepl.append(tof)
        velocities_wepl.append(velocity)

    return np.mean(wepls), np.mean(tofs_wepl), np.mean(velocities_wepl)


def tof_fit(
    wepls,
    output,
    verbose,
):

    results = []
    points = []

    with Pool() as pool:
        for wepl in wepls:
            result = pool.apply_async(
                process_wepl,
                (
                    wepl,
                    output,
                    verbose,
                ),
            )
            results.append(result)
        pool.close()
        pool.join()
    for result in results:
        point = result.get()
        points.append(point)

    np.savetxt(os.path.join(output, "data.txt"), points)

    data = np.array(points)
    wepls = data[:, 0]
    tofs = data[:, 1]
    velocities = data[:, 2]

    degs = range(3, 10)

    for deg in degs:
        tof_coeffs = np.polyfit(wepls, tofs, deg)
        vel_coeffs = np.polyfit(wepls, velocities, deg)
        np.savetxt(os.path.join(output, f"tof_coeffs_{deg}.txt"), tof_coeffs)
        np.savetxt(os.path.join(output, f"vel_coeffs_{deg}.txt"), vel_coeffs)


def process(args_info: argparse.Namespace):
    wepls = np.linspace(1.0, args_info.max_wepl, args_info.wepl_samples)

    results = []

    # Catch DeprecationWarnings coming from opengate and making PCT tests fail
    # To be removed once the warnings are gone from opengate
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)

        with Pool(maxtasksperchild=1) as pool:
            for wepl in wepls:
                result = pool.apply_async(
                    tof_fit_mc,
                    (
                        wepl,
                        args_info.output,
                        args_info.material,
                        args_info.phantom_width,
                        args_info.number_of_particles,
                        args_info.initial_energy,
                        args_info.seed,
                        args_info.verbose,
                    ),
                )
                results.append(result)
            pool.close()
            pool.join()
        for result in results:
            result.get()

    tof_fit(wepls, args_info.output, args_info.verbose)


def main(argv=None):
    parser = build_parser()
    args_info = parser.parse_args(argv)
    process(args_info)


if __name__ == "__main__":
    main()
