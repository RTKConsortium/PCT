#ifndef __pctProtonPairsIterator_h
#define __pctProtonPairsIterator_h

#include "itkImageRegionConstIterator.h"

namespace pct
{

template <typename TVectorType>
class ITK_TEMPLATE_EXPORT ProtonPairsIterator : public itk::ImageRegionConstIterator<itk::Image<TVectorType, 1>>
{
public:
  using ValueType = typename TVectorType::ValueType;
  using ImageType = itk::Image<TVectorType, 1>;
  using Superclass = itk::ImageRegionConstIterator<ImageType>;
  using RegionType = typename Superclass::RegionType;

  ProtonPairsIterator(const ImageType * img, const RegionType & region);

  itk::Vector<ValueType, 3>
  GetUpstreamPosition();

protected:
  int
  GetFieldIndex(const itk::MetaDataDictionary & metaDataDictionary, const std::string key);

  int m_upsteamPositionU;
  int m_upsteamPositionV;
  int m_upsteamPositionW;

}; // end of class

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#  include "pctProtonPairsIterator.hxx"
#endif

#endif
