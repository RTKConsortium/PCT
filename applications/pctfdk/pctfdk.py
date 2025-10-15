#!/usr/bin/env python
import argparse
import numpy as np
import itk
from itk import PCT as pct


def build_parser():
    parser = pct.PCTArgumentParser(
        description="Reconstruct a 3D volume from a sequence of projections [Feldkamp, David, Kress, 1984]."
    )
    # General
    parser.add_argument(
        "--verbose", "-v", help="Verbose execution", action="store_true"
    )
    parser.add_argument(
        "--geometry", "-g", help="XML geometry file name", type=str, required=True
    )
    parser.add_argument(
        "--path", "-p", help="Path containing projections", type=str, required=True
    )
    parser.add_argument(
        "--regexp",
        "-r",
        help="Regular expression to select projection files in path",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output", "-o", help="Output file name", type=str, required=True
    )
    parser.add_argument(
        "--lowmem",
        "-l",
        help="Load only one projection per thread in memory",
        action="store_true",
    )
    parser.add_argument(
        "--wpc",
        help="Water precorrection coefficients (default is no correction)",
        type=float,
        nargs="+",
    )

    # Ramp filter
    parser.add_argument(
        "--pad",
        help="Data padding parameter to correct for truncation",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--hann",
        help="Cut frequency for hann window in ]0,1] (0.0 disables it)",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--hannY",
        help="Cut frequency for hann window in ]0,1] (0.0 disables it)",
        type=float,
        default=0.0,
    )

    # Volume properties
    parser.add_argument(
        "--origin", help="Origin (default=centered)", type=float, nargs="+"
    )
    parser.add_argument(
        "--dimension", help="Deprecated. Use --size instead.", type=int, nargs="+"
    )
    parser.add_argument("--size", help="Size", type=int, nargs="+", default=[256])
    parser.add_argument(
        "--spacing", help="Spacing", type=float, nargs="+", default=[1.0]
    )
    parser.add_argument(
        "--direction", help="Direction (3x3 row-major)", type=float, nargs="+"
    )
    parser.add_argument(
        "--like",
        help="Copy info from image (origin, size, spacing, direction)",
        type=str,
    )
    return parser


def process(args_info: argparse.Namespace):
    from itk import RTK as rtk

    # Generate file names
    names = itk.RegularExpressionSeriesFileNames.New()
    names.SetDirectory(args_info.path)
    names.SetNumericSort(False)
    names.SetRegularExpression(args_info.regexp)
    names.SetSubMatch(0)
    if args_info.verbose:
        print(f"Regular expression matches {len(names.GetFileNames())} file(s)...")

    # Projections reader
    OutputPixelType = itk.F
    ProjectionImageType = itk.Image[OutputPixelType, 4]
    ReaderType = rtk.ProjectionsReader[ProjectionImageType]
    reader = ReaderType.New()
    reader.SetFileNames(names.GetFileNames())
    if args_info.wpc:
        reader.SetWaterPrecorrectionCoefficients([float(c) for c in args_info.wpc])

    # Geometry
    if args_info.verbose:
        print(f"Reading geometry information from {args_info.geometry}...")
    geometryReader = rtk.ThreeDCircularProjectionGeometryXMLFileReader.New()
    geometryReader.SetFilename(args_info.geometry)
    geometryReader.GenerateOutputInformation()

    # Short scan image filter
    pssf = pct.DDParkerShortScanImageFilter[ProjectionImageType].New()
    pssf.SetInput(reader.GetOutput())
    pssf.SetGeometry(geometryReader.GetOutputObject())
    pssf.InPlaceOff()

    # Create reconstructed image
    OutputImageType = itk.Image[OutputPixelType, 3]
    constantImageSource = rtk.ConstantImageSource[OutputImageType].New()
    rtk.SetConstantImageSourceFromArgParse(constantImageSource, args_info)

    # FDK reconstruction filtering
    FeldkampType = pct.FDKDDConeBeamReconstructionFilter[OutputImageType]
    feldkamp = FeldkampType.New()
    feldkamp.SetInput(0, constantImageSource.GetOutput())
    feldkamp.SetProjectionStack(pssf.GetOutput())
    feldkamp.SetGeometry(geometryReader.GetOutputObject())
    feldkamp.GetRampFilter().SetTruncationCorrection(args_info.pad)
    feldkamp.GetRampFilter().SetHannCutFrequency(args_info.hann)
    feldkamp.GetRampFilter().SetHannCutFrequencyY(args_info.hannY)

    # Write
    WriterType = itk.ImageFileWriter[OutputImageType]
    writer = WriterType.New()
    writer.SetFileName(args_info.output)
    writer.SetInput(feldkamp.GetOutput())

    if args_info.verbose:
        print("Reconstructing and writing... ", end="", flush=True)

    writerProbe = itk.TimeProbe()
    writerProbe.Start()
    writer.Update()
    writerProbe.Stop()

    if args_info.verbose:
        print(f"It took {writerProbe.GetMean()} {writerProbe.GetUnit()}")
        feldkamp.PrintTiming()


def main(argv=None):
    parser = build_parser()
    args_info = parser.parse_args(argv)
    process(args_info)


if __name__ == "__main__":
    main()
