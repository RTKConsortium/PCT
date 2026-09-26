#include "pctpocnewformat_ggo.h"
#include "pctProtonPairsIterator.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>

const std::string TEST_FILE = "/tmp/test.mhd";

int
main(int argc, char * argv[])
{
  GGO(pctpocnewformat, args_info); // RTK macro parsing options from .ggo file (rtkMacro.h)

  using VectorType = itk::Vector<float, 3>;
  using ImageType = itk::Image<VectorType, 1>;

  // Create a list-mode dataset that uses the new format

  auto                       img = ImageType::New();
  const ImageType::IndexType start = { { 0 } };
  const ImageType::SizeType  size = { { 10 } };
  ImageType::RegionType      region;
  region.SetSize(size);
  region.SetIndex(start);
  img->SetRegions(region);
  img->Allocate();

  itk::ImageRegionIterator<ImageType> it(img, img->GetLargestPossibleRegion());
  float                               p = 0.;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  { // Fill the data with dummy data
    float x = p + 1.;
    float y = p + 2.;
    float z = p + 3.;
    it.Set(itk::MakeVector(x, y, z));
    p++;
  }

  auto metaDataDictionary = img->GetMetaDataDictionary();
  itk::EncapsulateMetaData<int>(metaDataDictionary, "UpstreamPositionU", 0);
  itk::EncapsulateMetaData<int>(metaDataDictionary, "UpstreamPositionV", 1);
  itk::EncapsulateMetaData<int>(metaDataDictionary, "UpstreamPositionW", 2);
  img->SetMetaDataDictionary(metaDataDictionary);

  itk::WriteImage(img, TEST_FILE);

  // Read from a file that uses the new format

  auto img2 = itk::ReadImage<ImageType>(TEST_FILE);

  pct::ProtonPairsIterator<VectorType> it2(img2.GetPointer(), img2->GetLargestPossibleRegion());
  for (it2.GoToBegin(); !it2.IsAtEnd(); ++it2)
  {
    auto pos = it2.GetUpstreamPosition();
    std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
  }

  return EXIT_SUCCESS;
}
