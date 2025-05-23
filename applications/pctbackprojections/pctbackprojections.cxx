#include "pctbackprojections_ggo.h"
#include "rtkGgoFunctions.h"

#include "rtkThreeDCircularProjectionGeometryXMLFile.h"
#include "rtkProjectionsReader.h"
#include "pctFDKDDBackProjectionImageFilter.h"
#include "rtkJosephBackProjectionImageFilter.h"

#include <itkRegularExpressionSeriesFileNames.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbe.h>

int
main(int argc, char * argv[])
{
  GGO(pctbackprojections, args_info);

  using OutputPixelType = float;
  const unsigned int Dimension = 3;

  using ProjectionImageType = itk::Image<OutputPixelType, Dimension + 1>;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;

  // Geometry
  if (args_info.verbose_flag)
    std::cout << "Reading geometry information from " << args_info.geometry_arg << "..." << std::flush;
  rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReader;
  geometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(args_info.geometry_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(geometryReader->GenerateOutputInformation())
  if (args_info.verbose_flag)
    std::cout << " done." << std::endl;

  // Create an empty volume
  using ConstantImageSourceType = rtk::ConstantImageSource<OutputImageType>;
  ConstantImageSourceType::Pointer constantImageSource = ConstantImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctbackprojections>(constantImageSource,
                                                                                            args_info);

  // Generate file names
  itk::RegularExpressionSeriesFileNames::Pointer names = itk::RegularExpressionSeriesFileNames::New();
  names->SetDirectory(args_info.path_arg);
  names->SetNumericSort(false);
  names->SetRegularExpression(args_info.regexp_arg);
  names->SetSubMatch(0);

  if (args_info.verbose_flag)
    std::cout << "Reading " << names->GetFileNames().size() << " projection file(s)..." << std::flush;

  // Projections reader
  using ReaderType = rtk::ProjectionsReader<ProjectionImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames(names->GetFileNames());
  if (args_info.wpc_given)
  {
    std::vector<double> coeffs;
    coeffs.assign(args_info.wpc_arg, args_info.wpc_arg + args_info.wpc_given);
    reader->SetWaterPrecorrectionCoefficients(coeffs);
  }
  TRY_AND_EXIT_ON_ITK_EXCEPTION(reader->Update());

  if (args_info.verbose_flag)
    std::cout << " done." << std::endl;

  // Create back projection image filter
  if (args_info.verbose_flag)
    std::cout << "Backprojecting volume..." << std::flush;
  using BPType = pct::FDKDDBackProjectionImageFilter<OutputImageType, OutputImageType>;
  BPType::Pointer bp = BPType::New();
  bp->SetInput(constantImageSource->GetOutput());
  bp->SetProjectionStack(reader->GetOutput());
  bp->SetGeometry(geometryReader->GetOutputObject());
  TRY_AND_EXIT_ON_ITK_EXCEPTION(bp->Update())
  if (args_info.verbose_flag)
    std::cout << " done ." << std::endl;

  // Write
  if (args_info.verbose_flag)
    std::cout << "Writing... " << std::flush;
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(args_info.output_arg);
  writer->SetInput(bp->GetOutput());
  TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update());
  if (args_info.verbose_flag)
    std::cout << " done" << std::endl;

  return EXIT_SUCCESS;
}
