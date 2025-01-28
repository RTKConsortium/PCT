#include "pctProtonPairsToDistanceDrivenProjection.h"

#include "itkTestingMacros.h"

int
pctProtonPairsToDistanceDrivenProjectionTest(int argc, char* argv[])
{
  constexpr unsigned int Dimension = 3;
  using PixelType = float;
  using ImageType = itk::Image<PixelType, Dimension>;

  using FilterType = pct::ProtonPairsToDistanceDrivenProjection<ImageType, ImageType>;
  FilterType::Pointer filter = FilterType::New();

  ITK_EXERCISE_BASIC_OBJECT_METHODS(filter, ProtonPairsToDistanceDrivenProjection, InPlaceImageFilter);

  
  return EXIT_SUCCESS;
}
