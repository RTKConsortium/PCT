#ifndef __pctFDKDDConeBeamVarianceReconstructionFilter_h
#define __pctFDKDDConeBeamVarianceReconstructionFilter_h

#include "pctFFTVarianceImageFilter.h"
#include "pctFDKDDConeBeamReconstructionFilter.h"
#include "pctFFTVarianceImageFilter.h"

#include <itkExtractImageFilter.h>
#include <itkTimeProbe.h>

/** \class FDKDDConeBeamReconstructionFilter
 * TODO
 *
 * \author Jannis Dickmann
 */
namespace pct
{

template <class TInputImage, class TOutputImage = TInputImage, class TFFTPrecision = double>
class ITK_EXPORT FDKDDConeBeamVarianceReconstructionFilter
  : public pct::FDKDDConeBeamReconstructionFilter<TInputImage, TOutputImage>
{
public:
  typedef pct::FDKDDConeBeamReconstructionFilter<TInputImage, TOutputImage> Baseclass;

  /** Standard class typedefs. */
  typedef FDKDDConeBeamVarianceReconstructionFilter          Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                            Pointer;
  typedef itk::SmartPointer<const Self>                      ConstPointer;

  /** Some convenient typedefs. */
  typedef typename Baseclass::InputImageType  InputImageType;
  typedef typename Baseclass::OutputImageType OutputImageType;

  typedef typename Baseclass::ProjectionStackType    ProjectionStackType;
  typedef typename Baseclass::ProjectionStackPointer ProjectionStackPointer;

  /** Typedefs of each subfilter of this composite filter */
  typedef typename Baseclass::ExtractFilterType                                                ExtractFilterType;
  typedef typename Baseclass::WeightFilterType                                                 WeightFilterType;
  typedef pct::FFTVarianceImageFilter<ProjectionStackType, ProjectionStackType, TFFTPrecision> VarianceFilterType;
  typedef typename Baseclass::BackProjectionFilterType                                         BackProjectionFilterType;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkOverrideGetNameOfClassMacro(FDKDDConeBeamVarianceReconstructionFilter);

protected:
  FDKDDConeBeamVarianceReconstructionFilter();
  ~FDKDDConeBeamVarianceReconstructionFilter() {}

private:
  // purposely not implemented
  FDKDDConeBeamVarianceReconstructionFilter(const Self &);
  void
  operator=(const Self &);
}; // end of class

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#  include "pctFDKDDConeBeamVarianceReconstructionFilter.hxx"
#endif

#endif
