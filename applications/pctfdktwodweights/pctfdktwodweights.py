#!/usr/bin/env python
import argparse
import itk

from itk import PCT as pct


def build_parser():
    parser = pct.PCTArgumentParser(
        description="Weight projection series with 2D weights of the Feldkamp cone-beam reconstruction algorithm."
    )
    parser.add_argument(
        "--verbose", "-v", help="Verbose execution", action="store_true"
    )
    parser.add_argument(
        "--geometry", "-g", help="XML geometry file name", type=str, required=True
    )
    parser.add_argument(
        "--output", "-o", help="Output file name", type=str, required=True
    )
    parser.add_argument(
        "--divisions",
        "-d",
        help="Number of stream divisions to cope with large series",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--path", "-p", help="Path containing projections", required=True
    )
    parser.add_argument(
        "--regexp",
        "-r",
        help="Regular expression to select projection files in path",
        required=True,
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
    ProjectionImageType = itk.Image[itk.F, 4]
    ReaderType = rtk.ProjectionsReader[ProjectionImageType]
    reader = ReaderType.New()
    reader.SetFileNames(names.GetFileNames())
    reader.GenerateOutputInformation()

    # Geometry
    if args_info.verbose:
        print(f"Reading geometry information from {args_info.geometry}...")
    geometryReader = rtk.ThreeDCircularProjectionGeometryXMLFileReader.New()
    geometryReader.SetFilename(args_info.geometry)
    geometryReader.GenerateOutputInformation()

    # Weight filter
    WeightType = pct.FDKDDWeightProjectionFilter[ProjectionImageType]
    wf = WeightType.New()
    wf.SetInput(reader.GetOutput())
    wf.SetGeometry(geometryReader.GetOutputObject())
    wf.InPlaceOff()

    # Writer
    WriterType = itk.ImageFileWriter[ProjectionImageType]
    writer = WriterType.New()
    writer.SetFileName(args_info.output)
    writer.SetInput(wf.GetOutput())
    writer.SetNumberOfStreamDivisions(args_info.divisions)
    if args_info.verbose:
        print(f"Writing output to {args_info.output}...")
    writer.Update()


def main(argv=None):
    parser = build_parser()
    args_info = parser.parse_args(argv)
    process(args_info)


if __name__ == "__main__":
    main()
