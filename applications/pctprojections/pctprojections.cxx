#include "pctprojections_ggo.h"
#include "rtkMacro.h"

#include <rtkProjectionsReader.h>

#include <itkImageFileWriter.h>
#include <itkRegularExpressionSeriesFileNames.h>
#include <itkTimeProbe.h>

int
main(int argc, char * argv[])
{
  GGO(pctprojections, args_info);

  using OutputPixelType = float;
  const unsigned int Dimension = 4;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;

  // Generate file names
  itk::RegularExpressionSeriesFileNames::Pointer names = itk::RegularExpressionSeriesFileNames::New();
  names->SetDirectory(args_info.path_arg);
  names->SetNumericSort(false);
  names->SetRegularExpression(args_info.regexp_arg);
  names->SetSubMatch(0);

  // Projections reader
  using ReaderType = rtk::ProjectionsReader<OutputImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames(names->GetFileNames());
  if (args_info.wpc_given)
  {
    std::vector<double> coeffs;
    coeffs.assign(args_info.wpc_arg, args_info.wpc_arg + args_info.wpc_given);
    reader->SetWaterPrecorrectionCoefficients(coeffs);
  }

  // Write
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(args_info.output_arg);
  writer->SetInput(reader->GetOutput());
  TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->UpdateOutputInformation())
  writer->SetNumberOfStreamDivisions(1 + reader->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels() /
                                           (1024 * 1024 * 4));

  TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update())

  return EXIT_SUCCESS;
}
