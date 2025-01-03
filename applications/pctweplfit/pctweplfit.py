#!/usr/bin/env python
import argparse
import json
import sys
from multiprocessing import Pool, Manager

import matplotlib.pyplot as plt
import opengate as gate
import uproot
import numpy as np

DEFAULT_OUTPUT='pctweplfit'
DEFAULT_NUMBER_OF_PARTICLES=1e4
DEFAULT_PATH_TYPE='simple'
DEFAULT_PHANTOM_LENGTH_SAMPLES=10
DEFAULT_PHANTOM_WIDTH=40
DEFAULT_DETECTOR_DISTANCE=220.
DEFAULT_INITIAL_ENERGY=200.
DEFAULT_NUMBER_OF_DETECTORS=10
DEFAULT_POLYDEG_MIN=3
DEFAULT_POLYDEG_MAX=3
DEFAULT_VISU=False
DEFAULT_DISPLAY=False
DEFAULT_SAVEFIG=False
DEFAULT_VERBOSE=False

epsilon_mm = 1e-5

# Units
nm = gate.g4_units.nm
mm = gate.g4_units.mm
cm = gate.g4_units.cm
m = gate.g4_units.m
sec = gate.g4_units.second
MeV = gate.g4_units.MeV

def pv(verbose, *args, **kwargs):
    if verbose:
        print(*args, **kwargs)

def tof_fit_mc(
    phantom_length_mm,
    output=DEFAULT_OUTPUT,
    phantom_width_cm=DEFAULT_PHANTOM_WIDTH,
    detector_distance_mm=DEFAULT_DETECTOR_DISTANCE,
    number_of_particles=DEFAULT_NUMBER_OF_PARTICLES,
    initial_energy=DEFAULT_INITIAL_ENERGY,
    path_type=DEFAULT_PATH_TYPE,
    number_of_detectors=DEFAULT_NUMBER_OF_DETECTORS,
    visu=DEFAULT_VISU,
    verbose=DEFAULT_VERBOSE
):

    pv(verbose, "Starting simulation with following parameters: " + str(locals()))

    # Simulation
    sim = gate.Simulation()

    sim.random_engine = 'MersenneTwister'
    sim.random_seed = 'auto'
    sim.run_timing_intervals = [[0 * sec, 1 * sec]]
    sim.check_volumes_overlap = False
    sim.visu = visu
    sim.visu_type = 'vrml'
    sim.g4_verbose = False
    sim.progress_bar = verbose
    sim.number_of_threads = 1

    # Misc
    yellow = [1, 1, 0, 1]
    blue = [0, 0, 1, 1]

    # Geometry
    sim.world.material = 'G4_AIR'
    sim.world.size = [4 * m, 4 * m, 4 * m]
    sim.world.color = [0, 0, 0, 0]

    # Phantom
    if phantom_length_mm > 0.:
        phantom = sim.add_volume('Box', name='Phantom')
        phantom.size = [phantom_width_cm * cm, phantom_width_cm * cm, phantom_length_mm * mm]
        phantom.material = 'G4_WATER'
        phantom.color = blue
        phantom.set_max_step_size(1. * mm)

    # Beam
    source = sim.add_source('GenericSource', 'mybeam')
    source.particle = 'proton'
    source.energy.mono = initial_energy * MeV
    source.energy.type = 'mono'
    source.position.type = 'box'
    source.position.size = [1 * nm, 1 * nm, 1 * nm]
    source.position.translation = [0 * mm, 0 * mm, (-detector_distance_mm / 2 - 10) * mm]
    source.direction.type = 'momentum'
    source.direction.momentum = [0, 0, 1]
    source.n = number_of_particles

    # Physics list
    sim.physics_manager.physics_list_name = 'G4EmStandardPhysics_option4'

    # Phase spaces

    def add_detector(name, translation, attach_to_phantom=False):
        plane = sim.add_volume('Box', 'PlanePhaseSpace' + name)
        plane.size = [phantom_width_cm * cm, phantom_width_cm * cm, 1 * nm]
        plane.translation = translation
        plane.material = 'G4_AIR'
        plane.color = yellow
        if attach_to_phantom:
            plane.mother = phantom.name

        phase_space = sim.add_actor('PhaseSpaceActor', 'PhaseSpace' + name)
        phase_space.attached_to = plane.name
        phase_space.output_filename = f'{output}/l{int(phantom_length_mm)}_ps{name}.root'
        phase_space.attributes = [
            'EventID',
            'TrackID',
            'Position',
            'PreGlobalTime',
            'KineticEnergy'
        ]
        particle_filter = sim.add_filter('ParticleFilter', 'Filter' + name)
        particle_filter.particle = 'proton'

        phase_space.filters.append(particle_filter)

    add_detector('In', [0 * mm, 0 * mm, (-detector_distance_mm / 2 - epsilon_mm) * mm])
    add_detector('Out', [0 * mm, 0 * mm, (detector_distance_mm / 2 + epsilon_mm) * mm])

    if phantom_length_mm > 0. and path_type != 'phantom_length':
        for x in np.linspace(-phantom_length_mm / 2, phantom_length_mm / 2, number_of_detectors):
            add_detector(str(int(x)), [0 * mm, 0 * mm, x * mm], True)

    # Particle stats
    stat = sim.add_actor('SimulationStatisticsActor', 'stat')
    stat.output_filename = f'{output}/wepl{int(phantom_length_mm)}_stats.txt'

    sim.run()


def process_phantom_length(phantom_length, output, path_type, number_of_detectors, tofs, wepls, elosses, verbose):
    tofs_phantom_length = []
    wepls_phantom_length = []
    elosses_phantom_length = []

    data = uproot.concatenate(f'{output}/l{int(phantom_length)}_*.root', library='np')
    pv(verbose, "Loaded", len(data['EventID']), "events for phantom length", phantom_length)

    # Sort data if needed
    ws = data['Position_Z']
    if not np.all(ws[1:] > ws[:-1]):
        pv(verbose, "Sorting input data…")
        index_sorted = np.argsort(ws)
        for key in data.keys():
            data[key] = data[key][index_sorted]

    for n in np.unique(data['EventID']):
        event_mask = data['EventID'] == n

        number_of_hits = np.sum(event_mask)
        if number_of_hits == 0:
            continue
        if path_type != 'phantom_length' and number_of_hits < number_of_detectors + 2 and phantom_length > 0.:
            continue

        times = data['PreGlobalTime'][event_mask]
        tof = times[-1] - times[0]

        if phantom_length == 0.:
            wepl = 0.
        else:
            us = data['Position_X'][event_mask]
            vs = data['Position_Y'][event_mask]
            ws = data['Position_Z'][event_mask]
            if path_type == 'phantom_length':
                # Length of the phantom
                wepl = phantom_length
            elif path_type == 'simple':
                # Straight line between interaction position in first plane and in plane k
                wepl = np.sqrt((us[-2] - us[1])**2 + (vs[-2] - vs[1])**2 + (ws[-2] - ws[1])**2)
            elif path_type == 'realistic':
                # Path length through all detectors
                wepl = np.sum([
                        np.sqrt((us[w] - us[w - 1])**2 + (vs[w] - vs[w - 1])**2 + (ws[w] - ws[w - 1])**2)
                        for w in range(2, len(ws) - 1)
                    ])
            else:
                sys.exit(f"Invalid path time {path_type}!")

        eloss = data['KineticEnergy'][event_mask][0] - data['KineticEnergy'][event_mask][-1]

        tofs_phantom_length.append(tof)
        wepls_phantom_length.append(wepl)
        elosses_phantom_length.append(eloss)

    tofs[phantom_length] = tofs_phantom_length
    wepls[phantom_length] = wepls_phantom_length
    elosses[phantom_length] = elosses_phantom_length


def tof_fit(
    phantom_lengths,
    number_of_detectors,
    output=DEFAULT_OUTPUT,
    path_type=DEFAULT_PATH_TYPE,
    polydeg_min=DEFAULT_POLYDEG_MIN,
    polydeg_max=DEFAULT_POLYDEG_MAX,
    display=DEFAULT_DISPLAY,
    savefig=DEFAULT_SAVEFIG,
    verbose=DEFAULT_VERBOSE
):

    manager = Manager()
    tofs = manager.dict()
    wepls = manager.dict()
    elosses = manager.dict()

    results = []
    with Pool() as pool:
        for phantom_length in phantom_lengths:
            result = pool.apply_async(process_phantom_length, (phantom_length, output, path_type, number_of_detectors, tofs, wepls, elosses, verbose))
            results.append(result)
        pool.close()
        pool.join()
    for result in results:
        result.get()

    def fit(xs, ys, xlabel, ylabel):

        xmedians = [np.median(xs[phantom_length]) for phantom_length in phantom_lengths]
        ymedians = [np.median(ys[phantom_length]) for phantom_length in phantom_lengths]
        xpercentile25 = [np.percentile(xs[phantom_length], 25) for phantom_length in phantom_lengths]
        xpercentile75 = [np.percentile(xs[phantom_length], 75) for phantom_length in phantom_lengths]

        pv(verbose, f"Fitting {ylabel} to {xlabel}…")
        polydegs = range(polydeg_min, polydeg_max + 1)
        ps = {polydeg: np.polyfit(xmedians, ymedians, deg=polydeg).tolist() for polydeg in polydegs}
        pv(verbose, "Fitted coefficients:", ps)

        with open(f'{output}/{xlabel}_to_{ylabel}_medians.dat', 'w', encoding='utf-8') as f:
            for (xmedian, ymedian, xp25, xp75) in zip(xmedians, ymedians, xpercentile25, xpercentile75):
                f.write(f'{xmedian} {ymedian} {xp25} {xp75}\n')
        with open(f'{output}/{xlabel}_to_{ylabel}_points.dat', 'w', encoding='utf-8') as f:
            xpoints = [x for xx in xs.values() for x in xx]
            ypoints = [y for yy in ys.values() for y in yy]
            for (xpoint, ypoint) in zip(xpoints, ypoints):
                f.write(f'{xpoint} {ypoint}\n')

        if display or savefig:
            plt.figure()
            plt.plot(xmedians, ymedians, '+', label="Medians")

            tof_xs = np.linspace(np.min(xmedians), np.max(xmedians), 100)
            for d, p in ps.items():
                plt.plot(tof_xs, np.polyval(p, tof_xs), label=f"Polynomial fit (degree {d})")

            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.legend()

            if savefig:
                plt.savefig(f'{output}/{xlabel}_to_{ylabel}_fit.pdf')
            if display:
                plt.show()

        for d, p in ps.items():
            with open(f'{output}/{xlabel}_to_{ylabel}_fit_deg{d}.json', 'w', encoding='utf-8') as f:
                json.dump(p, f)

    fit(tofs, wepls, 'tof', 'wepl')
    fit(elosses, wepls, 'eloss', 'wepl')


def pctweplfit(
    output=DEFAULT_OUTPUT,
    number_of_particles=DEFAULT_NUMBER_OF_PARTICLES,
    path_type=DEFAULT_PATH_TYPE,
    phantom_length_samples=DEFAULT_PHANTOM_LENGTH_SAMPLES,
    phantom_width=DEFAULT_PHANTOM_WIDTH,
    detector_distance=DEFAULT_DETECTOR_DISTANCE,
    number_of_detectors=DEFAULT_NUMBER_OF_DETECTORS,
    initial_energy=DEFAULT_INITIAL_ENERGY,
    polydeg_min=DEFAULT_POLYDEG_MIN,
    polydeg_max=DEFAULT_POLYDEG_MAX,
    visu=DEFAULT_VISU,
    display=DEFAULT_DISPLAY,
    savefig=DEFAULT_SAVEFIG,
    verbose=DEFAULT_VERBOSE
):
    phantom_lengths = np.linspace(0., detector_distance, phantom_length_samples)

    results = []
    with Pool(maxtasksperchild=1) as pool:
        for phantom_length in phantom_lengths:
            result = pool.apply_async(tof_fit_mc, (
                phantom_length,
                output,
                phantom_width,
                detector_distance,
                number_of_particles,
                initial_energy,
                path_type,
                number_of_detectors,
                visu,
                verbose
            ))
            results.append(result)
        pool.close()
        pool.join()
    for result in results:
        result.get()

    tof_fit(phantom_lengths, number_of_detectors, output, path_type, polydeg_min, polydeg_max, display, savefig, verbose)


def main():

    parser = argparse.ArgumentParser(description="Convert TOF to WEPL using a fit on Monte Carlo data")
    parser.add_argument('-o', '--output', help="Path of outputs", default=DEFAULT_OUTPUT)
    parser.add_argument('-n', '--number-of-particles', help="Number of generated particles", default=DEFAULT_NUMBER_OF_PARTICLES, type=int)
    parser.add_argument('--path-type', help="How to compute proton path", choices=['phantom_length', 'simple', 'realistic'], default=DEFAULT_PATH_TYPE)
    parser.add_argument('--phantom-length-samples', help="Number of phantom length samples", default=DEFAULT_PHANTOM_LENGTH_SAMPLES, type=int)
    parser.add_argument('-w', '--phantom-width', help="Phantom width", default=DEFAULT_PHANTOM_WIDTH, type=float)
    parser.add_argument('-d', '--detector-distance', help="Distance between detectors", default=DEFAULT_DETECTOR_DISTANCE, type=float)
    parser.add_argument('-e', '--initial-energy', help="Initial energy of the protons (in MeV)", default=DEFAULT_INITIAL_ENERGY, type=float)
    parser.add_argument('--number-of-detectors', help="Number of detectors in the phantom", default=DEFAULT_NUMBER_OF_DETECTORS, type=int)
    parser.add_argument('--polydeg-min', help="Minimum polynom degree", default=DEFAULT_POLYDEG_MIN, type=int)
    parser.add_argument('--polydeg-max', help="Maximum polynom degree", default=DEFAULT_POLYDEG_MAX, type=int)
    parser.add_argument('--visu', help="Visualize Monte Carlo simulation", default=DEFAULT_VISU, action='store_true')
    parser.add_argument('--display', help="Display polynomial fit plot", default=DEFAULT_DISPLAY, action='store_true')
    parser.add_argument('--savefig', help="Write polynomial fit plot to disk", default=DEFAULT_SAVEFIG, action='store_true')
    parser.add_argument('--verbose', '-v', help="Verbose execution", default=DEFAULT_VERBOSE, action='store_true')
    args_info = parser.parse_args()

    pctweplfit(**vars(args_info))

if __name__ == '__main__':
    main()
