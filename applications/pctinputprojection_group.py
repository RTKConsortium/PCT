import itk
from itk import PCT as pct

__all__ = [
    "add_pctinputprojections_group",
    "GetProjectionsFileNamesFromArgParse",
]


# Mimicks pctinputprojections_section.ggo
def add_pctinputprojections_group(parser):
    pctinputprojections_group = parser.add_argument_group(
        "Input projections and their pre-processing"
    )
    pctinputprojections_group.add_argument(
        "--path", "-p", help="Path containing projections", required=True
    )
    pctinputprojections_group.add_argument(
        "--regexp",
        "-r",
        help="Regular expression to select projection files in path",
        required=True,
    )
    pctinputprojections_group.add_argument(
        "--nsort",
        help="Numeric sort for regular expression matches",
        action="store_true",
    )
    pctinputprojections_group.add_argument(
        "--submatch",
        help="Index of the submatch that will be used to sort matches",
        type=int,
        default=0,
    )
    pctinputprojections_group.add_argument(
        "--wpc",
        help="Water precorrection coefficients (default is no correction)",
        type=float,
        nargs="+",
    )


# Mimicks GetProjectionsFileNamesFromGgo
def GetProjectionsFileNamesFromArgParse(args_info):
    # Generate file names
    names = itk.RegularExpressionSeriesFileNames.New()
    names.SetDirectory(args_info.path)
    names.SetNumericSort(args_info.nsort)
    names.SetRegularExpression(args_info.regexp)
    names.SetSubMatch(args_info.submatch)

    if args_info.verbose:
        print(f"Regular expression matches {len(names.GetFileNames())} file(s)...")

    fileNames = []
    for fn in names.GetFileNames():
        imageio = itk.ImageIOFactory.CreateImageIO(
            fn, itk.CommonEnums.IOFileMode_ReadMode
        )
        if imageio is None:
            print(f"Ignoring file: {fn}")
            continue
        fileNames.append(fn)

    return fileNames


def SetProjectionsReaderFromArgParse(reader, args_info):
    fileNames = GetProjectionsFileNamesFromArgParse(args_info)

    # Water precorrection
    if args_info.wpc is not None:
        reader.SetWaterPrecorrectionCoefficients([float(c) for c in args_info.wpc])

    # Pass list to projections reader and update information
    reader.SetFileNames(fileNames)
    reader.UpdateOutputInformation()
