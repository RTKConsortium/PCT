#include "pctHoleFillingImageFilter.h"

#include <itkImageDuplicator.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkNeighborhoodIterator.h>
#include <itkNumericTraits.h>
#include <ostream>

namespace pct
{

template <typename TInputImage, typename TOutputImage>
void
HoleFillingImageFilter<TInputImage, TOutputImage>::GenerateData()
{
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using PixelType = typename OutputImageType::PixelType;

  // Allocate output as a copy of input
  auto initialDuplicator = itk::ImageDuplicator<OutputImageType>::New();
  initialDuplicator->SetInputImage(this->GetInput());
  initialDuplicator->Update();
  this->GetOutput()->Graft(initialDuplicator->GetOutput());

  unsigned int                       iteration = 0;
  typename OutputImageType::SizeType radius;
  radius.Fill(1);

  while (iteration < this->m_MaximumNumberOfIterations)
  {
    // Fast presence check: if there are no hole pixels left, exit early.
    bool                                           hasHole = false;
    itk::ImageRegionConstIterator<OutputImageType> checkIt(this->GetOutput(),
                                                           this->GetOutput()->GetLargestPossibleRegion());
    for (checkIt.GoToBegin(); !checkIt.IsAtEnd(); ++checkIt)
    {
      if (checkIt.Get() == this->m_HolePixelValue)
      {
        hasHole = true;
        break;
      }
    }
    if (!hasHole)
      break;

    auto duplicator = itk::ImageDuplicator<OutputImageType>::New();
    duplicator->SetInputImage(this->GetOutput());
    duplicator->Update();
    typename OutputImageType::Pointer snapshot = duplicator->GetOutput();

    itk::NeighborhoodIterator<OutputImageType> neighbor(radius, snapshot, snapshot->GetLargestPossibleRegion());
    itk::ImageRegionIterator<OutputImageType>  writer(this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());

    const unsigned int neighborSize = neighbor.Size();
    const unsigned int centerIndex = neighborSize / 2;

    for (neighbor.GoToBegin(), writer.GoToBegin(); !neighbor.IsAtEnd(); ++neighbor, ++writer)
    {
      PixelType centerPixel = neighbor.GetCenterPixel();
      if (centerPixel == this->m_HolePixelValue)
      {
        PixelType    sum = itk::NumericTraits<PixelType>::Zero;
        unsigned int valid = 0;

        for (unsigned int i = 0; i < neighborSize; ++i)
        {
          if (i == centerIndex)
            continue;
          bool      inBounds = false;
          PixelType p = neighbor.GetPixel(i, inBounds);
          if (inBounds && p != this->m_HolePixelValue)
          {
            sum += p;
            ++valid;
          }
        }

        if (valid > 0)
        {
          PixelType outValue = sum / valid;
          writer.Set(outValue);
        }
      }
    }
    ++iteration;
  }
}

} // namespace pct
