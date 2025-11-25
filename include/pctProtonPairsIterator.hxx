#include "pctProtonPairsIterator.h"

#include "itkMetaDataObject.h"

namespace pct
{

template <typename TVectorType>
ProtonPairsIterator<TVectorType>::ProtonPairsIterator(const ImageType * img, const RegionType & region)
  : itk::ImageRegionConstIterator<ImageType>(img, region)
{
  auto metaDataDictionary = img->GetMetaDataDictionary();
  this->m_upsteamPositionU = GetFieldIndex(metaDataDictionary, "UpstreamPositionU");
  this->m_upsteamPositionV = GetFieldIndex(metaDataDictionary, "UpstreamPositionV");
  this->m_upsteamPositionW = GetFieldIndex(metaDataDictionary, "UpstreamPositionW");
}

template <typename TVectorType>
int
ProtonPairsIterator<TVectorType>::GetFieldIndex(const itk::MetaDataDictionary & metaDataDictionary,
                                                const std::string               key)
{
  std::string metaData;
  if (itk::ExposeMetaData(metaDataDictionary, key, metaData))
  {
    // will throw a std::invalid_argument if number cannot be parsed
    return std::stoi(metaData);
  }
  return -1; // indicates that the key is not defined in the metadata dictionary
}

template <typename TVectorType>
itk::Vector<typename TVectorType::ValueType, 3>
ProtonPairsIterator<TVectorType>::GetUpstreamPosition()
{
  auto                      currentVec = this->Get();
  itk::Vector<ValueType, 3> vec;
  vec.SetNthComponent(0, currentVec[m_upsteamPositionU]);
  vec.SetNthComponent(1, currentVec[m_upsteamPositionV]);
  vec.SetNthComponent(2, currentVec[m_upsteamPositionW]);
  return vec;
}

} // namespace pct
