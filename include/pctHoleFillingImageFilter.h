#ifndef __pctHoleFillingImageFilter_h
#define __pctHoleFillingImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkSmartPointer.h>
#include <itkNumericTraits.h>

namespace pct
{
/**
 *
 * A simple iterative hole filling filter. A pixel equal to the HolePixel is
 * considered a hole. On each iteration, every hole pixel is replaced by the
 * average of its non-hole neighbors in a radius=1 neighborhood (center
 * excluded). Iterations repeat until no hole pixels remain or MaxIterations
 * is reached or an iteration makes no progress.
 *
 * This filter reproduces the behavior of the project's SmallHoleFiller
 * helper. For the original reference see:
 * https://www.insight-journal.org/browse/publication/835/
 */

template <typename TInputImage, typename TOutputImage = TInputImage>
class HoleFillingImageFilter : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  using Self = HoleFillingImageFilter;
  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);
  itkTypeMacro(HoleFillingImageFilter, ImageToImageFilter);

  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using PixelType = typename OutputImageType::PixelType;

  itkSetMacro(HolePixelValue, PixelType);
  itkGetConstMacro(HolePixelValue, PixelType);
  itkSetMacro(MaximumNumberOfIterations, unsigned int);
  itkGetConstMacro(MaximumNumberOfIterations, unsigned int);

protected:
  HoleFillingImageFilter() = default;
  ~HoleFillingImageFilter() override = default;

  void
  GenerateData() override;

private:
  ITK_DISALLOW_COPY_AND_MOVE(HoleFillingImageFilter);
  PixelType    m_HolePixelValue = itk::NumericTraits<PixelType>::Zero;
  unsigned int m_MaximumNumberOfIterations = 20;
};

} // namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#  include "pctHoleFillingImageFilter.hxx"
#endif

#endif
