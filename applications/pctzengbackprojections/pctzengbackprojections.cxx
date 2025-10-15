#include "pctzengbackprojections_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>

#include "pctZengWeightedBackProjectionImageFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int
main(int argc, char * argv[])
{
  GGO(pctzengbackprojections, args_info);

  using InputPixelType = float;
  const unsigned int Dimension = 4;
  using InputImageType = itk::Image<InputPixelType, Dimension>;

  auto reader = itk::ImageFileReader<InputImageType>::New();
  reader->SetFileName(args_info.input_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(reader->Update())

  using ZengFilterType = pct::ZengWeightedBackProjectionImageFilter<InputImageType>;
  ZengFilterType::Pointer zeng;
  zeng = ZengFilterType::New();
  zeng->SetInput(reader->GetOutput());
  TRY_AND_EXIT_ON_ITK_EXCEPTION(zeng->Update())

  auto writer = itk::ImageFileWriter<ZengFilterType::OutputImageType>::New();
  writer->SetInput(zeng->GetOutput(0));
  writer->SetFileName(args_info.outputc_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->UpdateLargestPossibleRegion())
  writer->SetInput(zeng->GetOutput(1));
  writer->SetFileName(args_info.outputs_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->UpdateLargestPossibleRegion())


  return EXIT_SUCCESS;
}
