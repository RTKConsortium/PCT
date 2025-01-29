#ifndef __pctDDParkerShortScanImageFilter_h
#define __pctDDParkerShortScanImageFilter_h

#include <itkInPlaceImageFilter.h>
#include <rtkThreeDCircularProjectionGeometry.h>
#include "rtkConfiguration.h"

/** \class DDParkerShortScanImageFilter
 *  \ingroup PCT
 * Weighting of image projections to handle off-centered panels
 * in tomography reconstruction. Based on [Wang, Med Phys, 2002].
 *
 * Note that the filter does nothing if the panel shift is less than 10%
 * of its size. Otherwise, it does the weighting described in the publication
 * and zero pads the data on the nearest side to the center.
 *
 * \author Simon Rit
 */
namespace pct
{

template <class TInputImage, class TOutputImage = TInputImage>
class ITK_EXPORT DDParkerShortScanImageFilter : public itk::InPlaceImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  using Self = DDParkerShortScanImageFilter;

  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;

  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Some convenient typedefs. */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using OutputImageRegionType = typename OutputImageType::RegionType;
  using WeightImageType = itk::Image<typename TOutputImage::PixelType, 1>;

  using GeometryType = rtk::ThreeDCircularProjectionGeometry;
  using GeometryPointer = GeometryType::Pointer;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkOverrideGetNameOfClassMacro(DDParkerShortScanImageFilter);

  /** Get / Set the object pointer to projection geometry */
  itkGetMacro(Geometry, GeometryPointer);
  itkSetMacro(Geometry, GeometryPointer);

protected:
  DDParkerShortScanImageFilter() { this->SetInPlace(true); }
  ~DDParkerShortScanImageFilter() {}

  virtual void
  DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread) override;

private:
  DDParkerShortScanImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  /** RTK geometry object */
  GeometryPointer m_Geometry;

  /** Superior and inferior position of the detector along the weighting direction, i.e. x.
   * The computed value account for the x projection offset of the geometry.
   */
  double m_InferiorCorner;
  double m_SuperiorCorner;
}; // end of class

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#  include "pctDDParkerShortScanImageFilter.hxx"
#endif

#endif
