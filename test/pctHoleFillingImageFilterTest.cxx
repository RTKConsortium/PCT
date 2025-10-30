#include "pctHoleFillingImageFilter.h"

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkTestingMacros.h"

int
pctHoleFillingImageFilterTest(int argc, char * argv[])
{
  using ImageType = itk::Image<int, 3>;

  auto                  img = ImageType::New();
  ImageType::RegionType region;
  region.SetIndex(itk::MakeIndex(0, 0, 0));
  region.SetSize(itk::MakeSize(8, 8, 8));

  img->SetRegions(region);
  img->Allocate();
  img->FillBuffer(42);

  // Insert some holes (use 0 as hole value)
  itk::ImageRegionIterator<ImageType> it(img, region);
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    ImageType::IndexType idx = it.GetIndex();
    if ((idx[0] + idx[1] + idx[2]) % 10 == 0)
    {
      it.Set(0);
    }
  }

  using FilterType = pct::HoleFillingImageFilter<ImageType>;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(img);
  filter->SetHolePixelValue(0);
  filter->SetMaximumNumberOfIterations(10);
  filter->Update();

  ImageType::Pointer out = filter->GetOutput();

  // Check no zeros remain and values are 42
  itk::ImageRegionIterator<ImageType> outIt(out, region);
  for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt)
  {
    ITK_TEST_EXPECT_EQUAL(outIt.Get(), 42);
  }

  std::cout << "\n\nTest PASSED! " << std::endl;

  return EXIT_SUCCESS;
}
