#include "pctfdk_ggo.h"
#include "rtkGgoFunctions.h"

#include "rtkThreeDCircularProjectionGeometryXMLFile.h"
#include "rtkProjectionsReader.h"
#include "pctDDParkerShortScanImageFilter.h"
#include "pctFDKDDConeBeamReconstructionFilter.h"

#include <itkRegularExpressionSeriesFileNames.h>
#include <itkImageFileWriter.h>

int
main(int argc, char * argv[])
{
  GGO(pctfdk, args_info);

  using OutputPixelType = float;
  const unsigned int Dimension = 3;

  using OutputImageType = itk::Image<OutputPixelType, Dimension>;

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(
    std::min<double>(8, itk::MultiThreaderBase::GetGlobalMaximumNumberOfThreads()));

  // Generate file names
  itk::RegularExpressionSeriesFileNames::Pointer names = itk::RegularExpressionSeriesFileNames::New();
  names->SetDirectory(args_info.path_arg);
  names->SetNumericSort(false);
  names->SetRegularExpression(args_info.regexp_arg);
  names->SetSubMatch(0);

  if (args_info.verbose_flag)
    std::cout << "Regular expression matches " << names->GetFileNames().size() << " file(s)..." << std::endl;

  // Projections reader
  using ProjectionImageType = itk::Image<OutputPixelType, Dimension + 1>;
  using ReaderType = rtk::ProjectionsReader<ProjectionImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames(names->GetFileNames());
  if (args_info.wpc_given)
  {
    std::vector<double> coeffs;
    coeffs.assign(args_info.wpc_arg, args_info.wpc_arg + args_info.wpc_given);
    reader->SetWaterPrecorrectionCoefficients(coeffs);
  }

  // Geometry
  if (args_info.verbose_flag)
    std::cout << "Reading geometry information from " << args_info.geometry_arg << "..." << std::endl;
  rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReader;
  geometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(args_info.geometry_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(geometryReader->GenerateOutputInformation())

  // Short scan image filter
  using PSSFType = pct::DDParkerShortScanImageFilter<ProjectionImageType>;
  PSSFType::Pointer pssf = PSSFType::New();
  pssf->SetInput(reader->GetOutput());
  pssf->SetGeometry(geometryReader->GetOutputObject());
  pssf->InPlaceOff();

  // Create reconstructed image
  using ConstantImageSourceType = rtk::ConstantImageSource<OutputImageType>;
  ConstantImageSourceType::Pointer constantImageSource = ConstantImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctfdk>(constantImageSource, args_info);

  // FDK reconstruction filtering
  using FDKCPUType = pct::FDKDDConeBeamReconstructionFilter<OutputImageType>;
  FDKCPUType::Pointer feldkamp = FDKCPUType::New();
  feldkamp->SetInput(0, constantImageSource->GetOutput());
  feldkamp->SetProjectionStack(pssf->GetOutput());
  feldkamp->SetGeometry(geometryReader->GetOutputObject());
  feldkamp->GetRampFilter()->SetTruncationCorrection(args_info.pad_arg);
  feldkamp->GetRampFilter()->SetHannCutFrequency(args_info.hann_arg);
  feldkamp->GetRampFilter()->SetHannCutFrequencyY(args_info.hannY_arg);

  // Write
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(args_info.output_arg);
  writer->SetInput(feldkamp->GetOutput());

  if (args_info.verbose_flag)
    std::cout << "Reconstructing and writing... " << std::flush;
  itk::TimeProbe writerProbe;

  writerProbe.Start();
  TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update());
  writerProbe.Stop();

  if (args_info.verbose_flag)
  {
    std::cout << "It took " << writerProbe.GetMean() << ' ' << writerProbe.GetUnit() << std::endl;
    feldkamp->PrintTiming(std::cout);
  }

  return EXIT_SUCCESS;
}
