#ifndef __pctFDKDDBackProjectionImageFilter_h
#define __pctFDKDDBackProjectionImageFilter_h

#include <itkInPlaceImageFilter.h>
#include <rtkBackProjectionImageFilter.h>

namespace pct
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT FDKDDBackProjectionImageFilter : public rtk::BackProjectionImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  using Self = FDKDDBackProjectionImageFilter;
  using Superclass = rtk::BackProjectionImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using ProjectionMatrixType = typename Superclass::ProjectionMatrixType;
  using OutputImageRegionType = typename TOutputImage::RegionType;
  using ProjectionImageType = TInputImage;
  using ProjectionImagePointer = typename ProjectionImageType::Pointer;
  using ProjectionPixelType = typename ProjectionImageType::PixelType;
  using ProjectionStackType = itk::Image<float, 4>;
  using ProjectionStackPointer = ProjectionStackType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(FDKDDBackProjectionImageFilter);

  virtual ProjectionImagePointer
  GetDDProjection(const unsigned int iProj);

  /** Get / Set the stack of projection images */
  itkGetMacro(ProjectionStack, ProjectionStackPointer);
  itkSetMacro(ProjectionStack, ProjectionStackPointer);

  /** Creates the #iProj index to index projection matrix with current inputs
      instead of the physical point to physical point projection matrix provided by Geometry */
  ProjectionMatrixType
  GetIndexToIndexProjectionMatrix(const unsigned int iProj, const ProjectionImageType * proj);

protected:
  FDKDDBackProjectionImageFilter()
  {
    this->SetNumberOfRequiredInputs(1);
    this->SetInPlace(true);
  };
  virtual ~FDKDDBackProjectionImageFilter() {};

  virtual void
  DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread) override;

  ProjectionStackPointer m_ProjectionStack;

private:
  FDKDDBackProjectionImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented
};

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#  include "pctFDKDDBackProjectionImageFilter.hxx"
#endif

#endif
