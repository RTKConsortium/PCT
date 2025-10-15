#ifndef __pctZengWeightedBackProjectionImageFilter_h
#define __pctZengWeightedBackProjectionImageFilter_h

#include <itkImageToImageFilter.h>
#include <rtkThreeDCircularProjectionGeometry.h>
#include "rtkConfiguration.h"

/** \class ZengWeightedBackProjectionImageFilter
 * \ingroup PCT
 * From an input 4D image where the 4th dimension is the angle, computes
 * the weighted backprojection for DBP described in [Zeng, Med Phys, 2007].
 *
 * \author Simon Rit
 */
namespace pct
{

template <class TInputImage,
          class TOutputImage = itk::Image<typename TInputImage::PixelType, TInputImage::ImageDimension - 1>>
class ITK_TEMPLATE_EXPORT ZengWeightedBackProjectionImageFilter
  : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  using Self = ZengWeightedBackProjectionImageFilter;
  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Some convenient typedefs. */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using OutputImageRegionType = typename TOutputImage::RegionType;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkOverrideGetNameOfClassMacro(ZengWeightedBackProjectionImageFilter);

protected:
  ZengWeightedBackProjectionImageFilter();
  ~ZengWeightedBackProjectionImageFilter() {}

  virtual void
  GenerateOutputInformation() override;
  virtual void
  DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread) override;

  using Superclass::MakeOutput;
  virtual itk::DataObject::Pointer
  MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType idx) override;

private:
  ZengWeightedBackProjectionImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented
}; // end of class

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#  include "pctZengWeightedBackProjectionImageFilter.hxx"
#endif

#endif
