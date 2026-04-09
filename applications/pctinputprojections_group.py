import itk
from itk import PCT as pct
import numpy as np

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
        "--newdirection",
        help="New value of input projections (before pre-processing)",
        type=float,
        nargs="+",
    )
    pctinputprojections_group.add_argument(
        "--neworigin",
        help="New origin of input projections (before pre-processing)",
        type=float,
        nargs="+",
    )
    pctinputprojections_group.add_argument(
        "--newspacing",
        help="New spacing of input projections (before pre-processing)",
        type=float,
        nargs="+",
    )
    pctinputprojections_group.add_argument(
        "--lowercrop",
        help="Lower boundary crop size",
        type=int,
        nargs="+",
        default=[0],
    )
    pctinputprojections_group.add_argument(
        "--uppercrop",
        help="Upper boundary crop size",
        type=int,
        nargs="+",
        default=[0],
    )
    pctinputprojections_group.add_argument(
        "--binning",
        help="Shrink / Binning factos in each direction",
        type=int,
        nargs="+",
        default=[1],
    )
    pctinputprojections_group.add_argument(
        "--wpc",
        help="Water precorrection coefficients (default is no correction)",
        type=float,
        nargs="+",
    )
    pctinputprojections_group.add_argument(
        "--radius",
        help="Radius of neighborhood for conditional median filtering",
        type=int,
        nargs="+",
        default=[0],
    )
    pctinputprojections_group.add_argument(
        "--multiplier",
        help="Threshold multiplier for conditional median filtering",
        type=float,
        default=0,
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

    # Vector component extraction (not in PCT ggo)

    # Change image information
    Dimension = reader.GetOutput().GetImageDimension()
    if args_info.newdirection is not None:
        direction = [args_info.newdirection[0]] * 9
        for i in range(min(9, len(args_info.newdirection))):
            direction[i] = args_info.newdirection[i]
        direction = np.array(direction).reshape((3, 3))
        reader.SetDirection(itk.matrix_from_array(direction))

    if args_info.newspacing is not None:
        spacing = itk.Vector[itk.D, Dimension]()
        spacing.Fill(args_info.newspacing[0])
        for i in range(len(args_info.newspacing)):
            spacing[i] = args_info.newspacing[i]
        reader.SetSpacing(spacing)

    if args_info.neworigin is not None:
        origin = itk.Point[itk.D, Dimension]()
        origin.Fill(args_info.neworigin[0])
        for i in range(len(args_info.neworigin)):
            origin[i] = args_info.neworigin[i]
        reader.SetOrigin(origin)

    # Crop boundaries
    upperCrop = [0] * Dimension
    lowerCrop = [0] * Dimension
    if args_info.lowercrop is not None:
        for i in range(len(args_info.lowercrop)):
            lowerCrop[i] = args_info.lowercrop[i]
    reader.SetLowerBoundaryCropSize(lowerCrop)
    if args_info.uppercrop is not None:
        for i in range(len(args_info.uppercrop)):
            upperCrop[i] = args_info.uppercrop[i]
    reader.SetUpperBoundaryCropSize(upperCrop)

    # Conditional median
    medianRadius = reader.GetMedianRadius()
    if args_info.radius is not None:
        for i in range(len(args_info.radius)):
            medianRadius[i] = args_info.radius[i]
    reader.SetMedianRadius(medianRadius)
    if args_info.multiplier is not None:
        reader.SetConditionalMedianThresholdMultiplier(args_info.multiplier)

    # Shrink / Binning
    binFactors = reader.GetShrinkFactors()
    if args_info.binning is not None:
        for i in range(len(args_info.binning)):
            binFactors[i] = args_info.binning[i]
    reader.SetShrinkFactors(binFactors)

    # Water precorrection
    if args_info.wpc is not None:
        reader.SetWaterPrecorrectionCoefficients(args_info.wpc)

    # Pass list to projections reader and update information
    reader.SetFileNames(fileNames)
    reader.UpdateOutputInformation()
