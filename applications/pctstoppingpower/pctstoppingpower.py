import argparse
import warnings
import itk
from itk import PCT as pct


def build_parser():
    parser = pct.PCTArgumentParser(
        description="Create a database of stopping power for a given material"
    )
    parser.add_argument("-o", "--output", help="Output file")
    parser.add_argument(
        "--max-energy", help="Maximum energy to simulate, in MeV", type=int, default=500
    )
    parser.add_argument(
        "-m", "--material", help="Material to consider", default="G4_WATER"
    )
    return parser


def process(args_info: argparse.Namespace):

    # Catch DeprecationWarnings coming from opengate and making PCT tests fail
    # To be removed once the warnings are gone from opengate
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)

        import opengate as gate

        cm = gate.g4_units.cm
        m = gate.g4_units.m

        sim = gate.Simulation()

        sim.world.material = "G4_AIR"
        sim.world.size = [6 * m, 6 * m, 6 * m]
        sim.world.color = [0, 0, 0, 0]

        samplebox = sim.add_volume("Box", "Samplebox")
        samplebox.size = [10 * cm, 10 * cm, 10 * cm]
        samplebox.translation = [0 * cm, 0 * cm, 2.5 * cm]
        samplebox.material = args_info.material

        sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"

        # EmCalc actor
        em_calc = sim.add_actor("EmCalculatorActor", "test")
        em_calc.attached_to = samplebox.name
        em_calc.is_ion = True
        em_calc.ion_params = "1 1"
        em_calc.material = args_info.material
        em_calc.nominal_energies = list(range(args_info.max_energy, 0, -1))
        em_calc.savefile_path = args_info.output

        sim.run()


def main(argv=None):
    parser = build_parser()
    args_info = parser.parse_args(argv)
    process(args_info)


if __name__ == "__main__":
    main()
