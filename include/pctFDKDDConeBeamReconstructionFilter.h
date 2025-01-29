#ifndef __pctFDKDDConeBeamReconstructionFilter_h
#define __pctFDKDDConeBeamReconstructionFilter_h

#include "pctFDKDDWeightProjectionFilter.h"
#include <rtkFFTRampImageFilter.h>
#include "pctFDKDDBackProjectionImageFilter.h"

#include <itkExtractImageFilter.h>
#include <itkTimeProbe.h>

/** \class FDKDDConeBeamReconstructionFilter
 * TODO
 *  \ingroup PCT
 * \author Simon Rit
 */
namespace pct
{

template <class TInputImage, class TOutputImage = TInputImage, class TFFTPrecision = double>
class ITK_EXPORT FDKDDConeBeamReconstructionFilter : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  using Self = FDKDDConeBeamReconstructionFilter;
  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Some convenient typedefs. */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;

  using ProjectionStackType = itk::Image<float, 4>;
  using ProjectionStackPointer = ProjectionStackType::Pointer;

  /** Typedefs of each subfilter of this composite filter */
  using ExtractFilterType = itk::ExtractImageFilter<ProjectionStackType, ProjectionStackType>;
  using WeightFilterType = pct::FDKDDWeightProjectionFilter<ProjectionStackType, ProjectionStackType>;
  using RampFilterType = rtk::FFTRampImageFilter<ProjectionStackType, ProjectionStackType, TFFTPrecision>;
  using BackProjectionFilterType = pct::FDKDDBackProjectionImageFilter<OutputImageType, OutputImageType>;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkOverrideGetNameOfClassMacro(FDKDDConeBeamReconstructionFilter);

  /** Get / Set the object pointer to projection geometry */
  virtual rtk::ThreeDCircularProjectionGeometry::Pointer
  GetGeometry();
  virtual void
  SetGeometry(const rtk::ThreeDCircularProjectionGeometry::Pointer _arg);

  /** Get pointer to the ramp filter used by the feldkamp reconstruction */
  typename RampFilterType::Pointer
  GetRampFilter()
  {
    return m_RampFilter;
  }

  void
  PrintTiming(std::ostream & os) const;

  /** Get / Set the stack of projection images */
  itkGetMacro(ProjectionStack, ProjectionStackPointer);
  itkSetMacro(ProjectionStack, ProjectionStackPointer);

protected:
  FDKDDConeBeamReconstructionFilter();
  ~FDKDDConeBeamReconstructionFilter() {}

  virtual void
  GenerateInputRequestedRegion() override;

  void
  GenerateOutputInformation() override;

  void
  GenerateData() override;

  /** The two inputs should not be in the same space so there is nothing
   * to verify. */
  virtual void
  VerifyInputInformation() const override
  {}

  /** Pointers to each subfilter of this composite filter */
  typename ExtractFilterType::Pointer        m_ExtractFilter;
  typename WeightFilterType::Pointer         m_WeightFilter;
  typename RampFilterType::Pointer           m_RampFilter;
  typename BackProjectionFilterType::Pointer m_BackProjectionFilter;

private:
  // purposely not implemented
  FDKDDConeBeamReconstructionFilter(const Self &);
  void
  operator=(const Self &);

  /** Probes to time reconstruction */
  itk::TimeProbe m_PreFilterProbe;
  itk::TimeProbe m_FilterProbe;
  itk::TimeProbe m_BackProjectionProbe;

  ProjectionStackPointer m_ProjectionStack;
}; // end of class

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#  include "pctFDKDDConeBeamReconstructionFilter.hxx"
#endif

#endif
