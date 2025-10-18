#!/usr/bin/env python
import argparse
import sys
from typing import Optional
import itk
from itk import PCT as pct


# Build a parser in the same style as pctpairprotons
def build_parser():
    parser = pct.PCTArgumentParser(
        description="Create the 3D sequence (2D + distance) of proton radiographies"
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Input file name containing the proton pairs",
        required=True,
    )
    parser.add_argument("-o", "--output", help="Output file name")
    parser.add_argument("--elosswepl", help="Output file name (alias for --output)")
    parser.add_argument(
        "-c", "--count", help="Image of count of proton pairs per pixel"
    )
    parser.add_argument(
        "--scatwepl", help="Image of scattering WEPL of proton pairs per pixel"
    )
    parser.add_argument("--noise", help="Image of WEPL variance per pixel")
    parser.add_argument(
        "-r",
        "--robust",
        help="Use robust estimation of scattering using 19.1 percentile.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-v", "--verbose", help="Verbose execution", action="store_true", default=False
    )

    parser.add_argument(
        "-s", "--source", help="Source position", type=float, default=0.0
    )
    parser.add_argument(
        "--quadricIn",
        help=(
            "Parameters of the entrance quadric support function, "
            "see http://education.siggraph.org/static/HyperGraph/raytrace/rtinter4.htm"
        ),
        type=float,
        nargs="+",
    )
    parser.add_argument(
        "--quadricOut",
        help=(
            "Parameters of the exit quadric support function, "
            "see http://education.siggraph.org/static/HyperGraph/raytrace/rtinter4.htm"
        ),
        type=float,
        nargs="+",
    )
    parser.add_argument(
        "--mlptype",
        help="Type of most likely path (schulte, polynomial, or krah)",
        choices=["schulte", "polynomial", "krah"],
        default="schulte",
    )
    parser.add_argument(
        "--mlptrackeruncert",
        help="Consider tracker uncertainties in MLP [Krah 2018, PMB]",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--mlppolydeg",
        help="Degree of the polynomial to approximate 1/beta^2p^2",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--ionpot",
        help="Ionization potential used in the reconstruction in eV",
        type=float,
        default=68.9984,
    )
    parser.add_argument(
        "--fill",
        help="Fill holes, i.e. pixels that were not hit by protons",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--trackerresolution", help="Tracker resolution in mm", type=float
    )
    parser.add_argument(
        "--trackerspacing", help="Tracker pair spacing in mm", type=float
    )
    parser.add_argument(
        "--materialbudget", help="Material budget x/X0 of tracker", type=float
    )

    parser.add_argument(
        "--origin", help="Origin (default=centered)", type=float, nargs="+"
    )
    parser.add_argument(
        "--dimension",
        help="Dimension(Deprecated) Use --size instead.",
        type=int,
        nargs="+",
        default=[256],
    )
    parser.add_argument("--size", help="Size", type=int, nargs="+", default=[256])
    parser.add_argument(
        "--spacing", help="Spacing", type=float, nargs="+", default=[1.0]
    )
    parser.add_argument("--direction", help="Direction", type=float, nargs="+")
    parser.add_argument(
        "--like",
        help="Copy information from this image (origin, size, spacing, direction)",
    )

    return parser


def process(args_info: argparse.Namespace):
    from itk import RTK as rtk

    if args_info.elosswepl and args_info.output:
        print("Only --output or --elosswepl should be provided", file=sys.stderr)
        sys.exit(1)

    OutputPixelType = itk.F
    Dimension = 3
    OutputImageType = itk.Image[OutputPixelType, Dimension]

    max_threads = int(min(8, itk.MultiThreaderBase.GetGlobalMaximumNumberOfThreads()))
    itk.MultiThreaderBase.SetGlobalMaximumNumberOfThreads(max_threads)

    # Create a stack of empty projection images
    constantImageSource = rtk.ConstantImageSource[OutputImageType].New()
    rtk.SetConstantImageSourceFromArgParse(constantImageSource, args_info)

    # Projection filter
    projection = pct.ProtonPairsToDistanceDrivenProjection[
        OutputImageType, OutputImageType
    ].New()
    projection.SetInput(constantImageSource.GetOutput())
    projection.SetProtonPairsFileName(args_info.input)
    projection.SetSourceDistance(args_info.source)
    projection.SetMostLikelyPathType(args_info.mlptype)
    projection.SetMostLikelyPathPolynomialDegree(args_info.mlppolydeg)
    projection.SetMostLikelyPathTrackerUncertainties(args_info.mlptrackeruncert)
    if args_info.trackerresolution is not None:
        projection.SetTrackerResolution(args_info.trackerresolution)
    if args_info.trackerspacing is not None:
        projection.SetTrackerPairSpacing(args_info.trackerspacing)
    if args_info.materialbudget is not None:
        projection.SetMaterialBudget(args_info.materialbudget)
    projection.SetIonizationPotential(args_info.ionpot * 1e-6)
    projection.SetRobust(args_info.robust)
    projection.SetComputeScattering(bool(args_info.scatwepl))
    projection.SetComputeNoise(bool(args_info.noise))

    if args_info.quadricIn:
        # quadric = object surface
        qIn = rtk.QuadricShape.New()
        qIn.SetA(args_info.quadricIn[0])
        qIn.SetB(args_info.quadricIn[1])
        qIn.SetC(args_info.quadricIn[2])
        qIn.SetD(args_info.quadricIn[3])
        qIn.SetE(args_info.quadricIn[4])
        qIn.SetF(args_info.quadricIn[5])
        qIn.SetG(args_info.quadricIn[6])
        qIn.SetH(args_info.quadricIn[7])
        qIn.SetI(args_info.quadricIn[8])
        qIn.SetJ(args_info.quadricIn[9])
        projection.SetQuadricIn(qIn)
    if args_info.quadricOut:
        qOut = rtk.QuadricShape.New()
        qOut.SetA(args_info.quadricOut[0])
        qOut.SetB(args_info.quadricOut[1])
        qOut.SetC(args_info.quadricOut[2])
        qOut.SetD(args_info.quadricOut[3])
        qOut.SetE(args_info.quadricOut[4])
        qOut.SetF(args_info.quadricOut[5])
        qOut.SetG(args_info.quadricOut[6])
        qOut.SetH(args_info.quadricOut[7])
        qOut.SetI(args_info.quadricOut[8])
        qOut.SetJ(args_info.quadricOut[9])
        projection.SetQuadricOut(qOut)

    projection.Update()

    filler = pct.SmallHoleFiller[OutputImageType]()
    if args_info.fill:
        filler.SetImage(projection.GetOutput())
        filler.SetHolePixel(0.0)
        filler.Fill()

    cii = itk.ChangeInformationImageFilter[OutputImageType].New()
    if args_info.fill:
        cii.SetInput(filler.GetOutput())
    else:
        cii.SetInput(projection.GetOutput())
    cii.ChangeOriginOn()
    cii.ChangeDirectionOn()
    cii.ChangeSpacingOn()
    cii.SetOutputDirection(projection.GetOutput().GetDirection())
    cii.SetOutputOrigin(projection.GetOutput().GetOrigin())
    cii.SetOutputSpacing(projection.GetOutput().GetSpacing())

    if args_info.elosswepl or args_info.output:
        # Write
        out_name = args_info.elosswepl or args_info.output
        itk.imwrite(cii.GetOutput(), out_name)

    if args_info.count:
        # Write
        itk.imwrite(projection.GetCount(), args_info.count)

    if args_info.scatwepl:
        # Write
        itk.imwrite(projection.GetAngle(), args_info.scatwepl)

    if args_info.noise:
        # Write
        itk.imwrite(projection.GetSquaredOutput(), args_info.noise)


def main(argv: Optional[list[str]] = None):
    parser = build_parser()
    args_info = parser.parse_args(argv)
    process(args_info)


if __name__ == "__main__":
    main()
