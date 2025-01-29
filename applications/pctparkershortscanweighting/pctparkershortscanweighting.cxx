#include "pctparkershortscanweighting_ggo.h"
#include "rtkGgoFunctions.h"

#include <rtkThreeDCircularProjectionGeometryXMLFile.h>
#include <rtkProjectionsReader.h>
#include "pctDDParkerShortScanImageFilter.h"

#include <itkImageFileWriter.h>
#include <itkRegularExpressionSeriesFileNames.h>

int
main(int argc, char * argv[])
{
  GGO(pctparkershortscanweighting, args_info);

  using OutputPixelType = float;
  const unsigned int Dimension = 4;

  using OutputImageType = itk::Image<OutputPixelType, Dimension>;

  // Generate file names
  itk::RegularExpressionSeriesFileNames::Pointer names = itk::RegularExpressionSeriesFileNames::New();
  names->SetDirectory(args_info.path_arg);
  names->SetNumericSort(false);
  names->SetRegularExpression(args_info.regexp_arg);
  names->SetSubMatch(0);

  if (args_info.verbose_flag)
    std::cout << "Regular expression matches " << names->GetFileNames().size() << " file(s)..." << std::endl;

  // Projections reader
  using ReaderType = rtk::ProjectionsReader<OutputImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames(names->GetFileNames());
  TRY_AND_EXIT_ON_ITK_EXCEPTION(reader->GenerateOutputInformation())

  // Geometry
  if (args_info.verbose_flag)
    std::cout << "Reading geometry information from " << args_info.geometry_arg << "..." << std::endl;
  rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReader;
  geometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(args_info.geometry_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(geometryReader->GenerateOutputInformation())

  // Short scan image filter
  using PSSFType = pct::DDParkerShortScanImageFilter<OutputImageType>;
  PSSFType::Pointer pssf = PSSFType::New();
  pssf->SetInput(reader->GetOutput());
  pssf->SetGeometry(geometryReader->GetOutputObject());
  pssf->InPlaceOff();

  // Write
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(args_info.output_arg);
  writer->SetInput(pssf->GetOutput());
  writer->SetNumberOfStreamDivisions(args_info.divisions_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update())

  return EXIT_SUCCESS;
}
