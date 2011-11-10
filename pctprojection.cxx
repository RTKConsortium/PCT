#include "pctprojection_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>
#include <itkConstantImageSource.h>

#include "itkProtonPairsToDistanceDrivenProjection.h"

#include <itkImageFileWriter.h>
#include <itkRegularExpressionSeriesFileNames.h>
#include <itkTimeProbe.h>

int main(int argc, char * argv[])
{
  GGO(pctprojection, args_info);

  typedef float OutputPixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  // Create a stack of empty projection images
  typedef itk::ConstantImageSource< OutputImageType > ConstantImageSourceType;
  ConstantImageSourceType::Pointer constantImageSource = ConstantImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctprojection>(constantImageSource, args_info);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( constantImageSource->Update() );

  // Projection filter
  typedef itk::ProtonPairsToDistanceDrivenProjection<OutputImageType, OutputImageType> ProjectionFilter;
  ProjectionFilter::Pointer projection = ProjectionFilter::New();
  projection->SetInput( constantImageSource->GetOutput() );
  projection->SetProtonPairsFileName( args_info.input_arg );
  projection->SetSourceDistance( args_info.source_arg );

  if(args_info.quadricIn_given)
    {
    ProjectionFilter::RQIType::Pointer qIn = ProjectionFilter::RQIType::New();
    qIn->SetA(args_info.quadricIn_arg[0]);
    qIn->SetB(args_info.quadricIn_arg[1]);
    qIn->SetC(args_info.quadricIn_arg[2]);
    qIn->SetD(args_info.quadricIn_arg[3]);
    qIn->SetE(args_info.quadricIn_arg[4]);
    qIn->SetF(args_info.quadricIn_arg[5]);
    qIn->SetG(args_info.quadricIn_arg[6]);
    qIn->SetH(args_info.quadricIn_arg[7]);
    qIn->SetI(args_info.quadricIn_arg[8]);
    qIn->SetJ(args_info.quadricIn_arg[9]);
    projection->SetQuadricIn(qIn);
    }
  if(args_info.quadricOut_given)
    {
    ProjectionFilter::RQIType::Pointer qOut = ProjectionFilter::RQIType::New();
    qOut->SetA(args_info.quadricOut_arg[0]);
    qOut->SetB(args_info.quadricOut_arg[1]);
    qOut->SetC(args_info.quadricOut_arg[2]);
    qOut->SetD(args_info.quadricOut_arg[3]);
    qOut->SetE(args_info.quadricOut_arg[4]);
    qOut->SetF(args_info.quadricOut_arg[5]);
    qOut->SetG(args_info.quadricOut_arg[6]);
    qOut->SetH(args_info.quadricOut_arg[7]);
    qOut->SetI(args_info.quadricOut_arg[8]);
    qOut->SetJ(args_info.quadricOut_arg[9]);
    projection->SetQuadricOut(qOut);
    }

  TRY_AND_EXIT_ON_ITK_EXCEPTION( projection->Update() );

  // Write
  typedef itk::ImageFileWriter<  OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( args_info.output_arg );
  writer->SetInput( projection->GetOutput() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() )

  return EXIT_SUCCESS;
}
