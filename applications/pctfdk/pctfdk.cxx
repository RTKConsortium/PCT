#include "pctfdk_ggo.h"
#include "pctGgoFunctions.h"
#include "rtkGgoFunctions.h"

#include "rtkThreeDCircularProjectionGeometryXMLFile.h"
#include "rtkProjectionsReader.h"
#include "pctDDParkerShortScanImageFilter.h"
#include "pctFDKDDConeBeamReconstructionFilter.h"
#include "pctFDKDDConeBeamVarianceReconstructionFilter.h"

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

  // Projections reader
  using ProjectionImageType = itk::Image<OutputPixelType, Dimension + 1>;
  auto reader = rtk::ProjectionsReader<ProjectionImageType>::New();
  pct::SetProjectionsReaderFromGgo<rtk::ProjectionsReader<ProjectionImageType>, args_info_pctfdk>(reader, args_info);

  // Geometry
  if (args_info.verbose_flag)
    std::cout << "Reading geometry information from " << args_info.geometry_arg << "..." << std::endl;
  auto geometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(args_info.geometry_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(geometryReader->GenerateOutputInformation())

  // Short scan image filter
  auto pssf = pct::DDParkerShortScanImageFilter<ProjectionImageType>::New();
  pssf->SetInput(reader->GetOutput());
  pssf->SetGeometry(geometryReader->GetOutputObject());
  pssf->InPlaceOff();

  // Create reconstructed image
  using ConstantImageSourceType = rtk::ConstantImageSource<OutputImageType>;
  auto constantImageSource = ConstantImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctfdk>(constantImageSource, args_info);

  // FDK reconstruction filtering
  using FDKCPUType = pct::FDKDDConeBeamReconstructionFilter<OutputImageType>;
  FDKCPUType::Pointer feldkamp = nullptr;
  if (args_info.variance_flag)
  {
    feldkamp = pct::FDKDDConeBeamVarianceReconstructionFilter<OutputImageType>::New();
  }
  else
  {
    feldkamp = FDKCPUType::New();
  }

  feldkamp->SetInput(0, constantImageSource->GetOutput());
  feldkamp->SetProjectionStack(pssf->GetOutput());
  feldkamp->SetGeometry(geometryReader->GetOutputObject());
  feldkamp->GetRampFilter()->SetTruncationCorrection(args_info.pad_arg);
  feldkamp->GetRampFilter()->SetHannCutFrequency(args_info.hann_arg);
  feldkamp->GetRampFilter()->SetHannCutFrequencyY(args_info.hannY_arg);

  // Write
  auto writer = itk::ImageFileWriter<OutputImageType>::New();
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
    feldkamp->PrintTiming();
  }

  return EXIT_SUCCESS;
}
