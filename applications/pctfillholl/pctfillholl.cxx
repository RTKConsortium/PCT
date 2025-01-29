#include "pctfillholl_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>

#include "SmallHoleFiller.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkChangeInformationImageFilter.h>

int
main(int argc, char * argv[])
{
  GGO(pctfillholl, args_info);

  using OutputPixelType = float;
  const unsigned int Dimension = 3;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;

  // Reader
  using ReaderType = itk::ImageFileReader<OutputImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(args_info.input_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(reader->Update());

  SmallHoleFiller<OutputImageType> filler;
  filler.SetImage(reader->GetOutput());
  filler.SetHolePixel(0.);
  filler.Fill();

  using CIIType = itk::ChangeInformationImageFilter<OutputImageType>;
  CIIType::Pointer cii = CIIType::New();
  cii->SetInput(filler.GetOutput());
  cii->ChangeOriginOn();
  cii->ChangeDirectionOn();
  cii->ChangeSpacingOn();
  cii->SetOutputDirection(reader->GetOutput()->GetDirection());
  cii->SetOutputOrigin(reader->GetOutput()->GetOrigin());
  cii->SetOutputSpacing(reader->GetOutput()->GetSpacing());

  // Write
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(args_info.output_arg);
  writer->SetInput(cii->GetOutput());
  TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update());

  return EXIT_SUCCESS;
}
