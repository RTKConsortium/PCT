import argparse

from itk import PCT as pct
from itk import RTK as rtk
import itk

def main():
  
  # Parse arguments
  parser = argparse.ArgumentParser(description="Weight projection serie with 2D weights of Feldkamp cone-beam reconstruction algorithm.")
  parser.add_argument("--geometry", "-g", help="XML geometry file name", required=True)
  parser.add_argument("--path", "-p", help="Path containing projections", required=True)
  parser.add_argument("--regexp", "-r", help="Regular expression to select projection files in path", required=True)
  parser.add_argument("--output", "-o", help="Output file name", required=True)
  parser.add_argument("--divisions", "-d", help="Number of stream divisions to cope with large serie", required=False, default=1, type=int)
  parser.add_argument("--verbose", "-v", help="Verbose execution", type=bool)
  args=parser.parse_args()
  
  # Generate file names
  names = itk.RegularExpressionSeriesFileNames.New()
  names.SetDirectory(args.path)
  names.SetNumericSort(False)
  names.SetRegularExpression(args.regexp)
  names.SetSubMatch(0)

  if args.verbose: print(f"Regular expression matches {len(names.GetFileNames())} file(s)...")

  # Projections reader
  ProjectionImageType = itk.Image[itk.F, 4]
  ReaderType = rtk.ProjectionsReader[ProjectionImageType]
  reader = ReaderType.New()
  reader.SetFileNames(names.GetFileNames())
  reader.GenerateOutputInformation()

  # Geometry
  if args.verbose: print(f"Reading geometry information from {args.geometry}...")
  geometryReader = rtk.ThreeDCircularProjectionGeometryXMLFileReader.New()
  geometryReader.SetFilename(args.geometry)
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
  writer.SetFileName(args.output)
  writer.SetInput(wf.GetOutput())
  writer.SetNumberOfStreamDivisions(args.divisions)
  writer.Update()


if __name__ == '__main__':
  main()