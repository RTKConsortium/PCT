/*=========================================================================
 *
 *  Copyright RTK Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef pctFFTVarianceImageFilter_hxx
#define pctFFTVarianceImageFilter_hxx

namespace pct
{

template <class TInputImage, class TOutputImage, class TFFTPrecision>
FFTVarianceImageFilter<TInputImage, TOutputImage, TFFTPrecision>::FFTVarianceImageFilter()
{
  //  this->m_HannCutFrequency = 0.;
  //  this->m_CosineCutFrequency = 0.;
  //  this->m_HammingFrequency = 0.;
  //  this->m_HannCutFrequency = 0.;
  //  this->m_HannCutFrequencyY = 0.;
  //  this->m_RamLakCutFrequency = 0.;
  //  this->m_SheppLoganCutFrequency = 0.;
}

template <class TInputImage, class TOutputImage, class TFFTPrecision>
void
FFTVarianceImageFilter<TInputImage, TOutputImage, TFFTPrecision>::UpdateFFTProjectionsConvolutionKernel(
  const SizeType s)
{
  // if(this->m_KernelFFT.GetPointer() != ITK_NULLPTR && s == this->m_PreviousKernelUpdateSize)
  //   {
  //   return;
  //   }
  // this->m_PreviousKernelUpdateSize = s;

  const int width = s[0];
  const int height = s[1];

  // Allocate kernel
  SizeType size;
  size.Fill(1);
  size[0] = width;
  FFTInputImagePointer kernel = FFTInputImageType::New();
  kernel->SetRegions(size);
  kernel->Allocate();
  kernel->FillBuffer(0.);

  // Compute kernel in space domain (see Kak & Slaney, chapter 3 equation 61
  // page 72) although spacing is not squared according to equation 69 page 75
  double    spacing = this->GetInput()->GetSpacing()[0];
  IndexType ix, jx;
  ix.Fill(0);
  jx.Fill(0);
  kernel->SetPixel(ix, 1. / (4. * spacing));
  for (ix[0] = 1, jx[0] = size[0] - 1; ix[0] < typename IndexType::IndexValueType(size[0] / 2); ix[0] += 2, jx[0] -= 2)
  {
    double v = ix[0] * vnl_math::pi;
    v = -1. / (v * v * spacing);
    kernel->SetPixel(ix, v);
    kernel->SetPixel(jx, v);
  }

  //  /////////////////////////////////////////////////////////////////////////////////////////////
  //  // DEBUG
  //  typedef itk::ImageRegionIteratorWithIndex<FFTInputImageType> InverseIteratorType2;
  //  InverseIteratorType2 itiK2(kernel, kernel->GetLargestPossibleRegion());
  //  std::ofstream out0("kernel_real.csv");
  //  itiK2.GoToBegin();
  //  for(; !itiK2.IsAtEnd(); ++itiK2)
  //  {
  //    out0 << itiK2.Get() << std::endl;
  //  }
  //  out0.close();
  //  // DEBUG

  // FFT kernel
  typedef itk::RealToHalfHermitianForwardFFTImageFilter<FFTInputImageType, FFTOutputImageType> FFTType;
  typename FFTType::Pointer                                                                    fftK = FFTType::New();
  fftK->SetInput(kernel);
  fftK->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
  fftK->Update();
  this->m_KernelFFT = fftK->GetOutput();

  // Windowing (if enabled)
  typedef itk::ImageRegionIteratorWithIndex<typename FFTType::OutputImageType> IteratorType;
  IteratorType itK(this->m_KernelFFT, this->m_KernelFFT->GetLargestPossibleRegion());

  unsigned int n = this->m_KernelFFT->GetLargestPossibleRegion().GetSize(0);

  itK.GoToBegin();
  if (this->GetHannCutFrequency() > 0.)
  {
    const unsigned int ncut = itk::Math::Round<double>(n * std::min(1.0, this->GetHannCutFrequency()));
    for (unsigned int i = 0; i < ncut; i++, ++itK)
    {
      double k = 0.5 * (1 + std::cos(itk::Math::pi * i / ncut));
      itK.Set(itK.Get() * TFFTPrecision(k));
    }
  }
  else if (this->GetCosineCutFrequency() > 0.)
  {
    const unsigned int ncut = itk::Math::Round<double>(n * std::min(1.0, this->GetCosineCutFrequency()));
    for (unsigned int i = 0; i < ncut; i++, ++itK)
    {
      double k = std::cos(0.5 * itk::Math::pi * i / ncut);
      itK.Set(itK.Get() * TFFTPrecision(k));
    }
  }
  else if (this->GetHammingFrequency() > 0.)
  {
    const unsigned int ncut = itk::Math::Round<double>(n * std::min(1.0, this->GetHammingFrequency()));
    for (unsigned int i = 0; i < ncut; i++, ++itK)
    {
      double k = 0.54 + 0.46 * (std::cos(itk::Math::pi * i / ncut));
      itK.Set(itK.Get() * TFFTPrecision(k));
    }
  }
  else if (this->GetRamLakCutFrequency() > 0.)
  {
    const unsigned int ncut = itk::Math::Round<double>(n * std::min(1.0, this->GetRamLakCutFrequency()));
    for (unsigned int i = 0; i < ncut; i++, ++itK)
    {
    }
  }
  else if (this->GetSheppLoganCutFrequency() > 0.)
  {
    const unsigned int ncut = itk::Math::Round<double>(n * std::min(1.0, this->GetSheppLoganCutFrequency()));
    // sinc(0) --> is 1
    ++itK;
    for (unsigned int i = 1; i < ncut; i++, ++itK)
    {
      double x = 0.5 * vnl_math::pi * i / ncut;
      double k = std::sin(x) / x;
      itK.Set(itK.Get() * TFFTPrecision(k));
    }
  }
  else
  {
    itK.GoToReverseBegin();
    ++itK;
  }

  for (; !itK.IsAtEnd(); ++itK)
  {
    itK.Set(itK.Get() * TFFTPrecision(0.));
  }

  //  /////////////////////////////////////////////////////////////////////////////////////////////
  //  // DEBUG
  //  std::ofstream out1("kernel_fourier.csv");
  //  itK.GoToBegin();
  //  for(; !itK.IsAtEnd(); ++itK)
  //  {
  //    out1 << itK.Get().real() << " " << itK.Get().imag() << std::endl;
  //  }
  //  out1.close();
  //  // DEBUG

  // iFFT kernel
  typedef itk::HalfHermitianToRealInverseFFTImageFilter<FFTOutputImageType, FFTInputImageType> IFFTType;
  typename IFFTType::Pointer                                                                   ifftK = IFFTType::New();
  ifftK->SetInput(this->m_KernelFFT);
  ifftK->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
  ifftK->Update();
  typename FFTType::InputImageType::Pointer KernelIFFT = ifftK->GetOutput();

  // calculate ratio g_C/g^2
  typedef itk::ImageRegionIteratorWithIndex<typename FFTType::InputImageType> InverseIteratorType;
  InverseIteratorType                         itiK(KernelIFFT, KernelIFFT->GetLargestPossibleRegion());
  typename FFTType::InputImageType::PixelType sumgc = 0.;
  typename FFTType::InputImageType::PixelType sumg2 = 0.;
  typename FFTType::InputImageType::IndexType idxshifted;
  const unsigned int                          maxidx = KernelIFFT->GetLargestPossibleRegion().GetSize()[0];
  for (; !itiK.IsAtEnd(); ++itiK)
  {
    idxshifted = itiK.GetIndex();
    if (idxshifted[0] == 0)
      idxshifted[0] = maxidx - 1;
    else
      idxshifted[0] -= 1;

    sumgc += itiK.Get() * KernelIFFT->GetPixel(idxshifted);
    sumg2 += itiK.Get() * itiK.Get();
  }
  const typename FFTType::InputImageType::PixelType ratiogcg2 = sumgc / sumg2;

  // numerical integration to calculate f_interp
  const double aprecision = 0.00001;
  double       finterp = 0.;
  for (unsigned int i = 0; i < int(1. / aprecision); i++)
  {
    const double a = double(i) * aprecision;
    finterp += (1 - a) * (1 - a) + 2 * ratiogcg2 * (1 - a) * a + a * a;
  }
  finterp *= aprecision;

  //  /////////////////////////////////////////////////////////////////////////////////////////////
  //  // DEBUG
  //  std::ofstream out9("kernel_window_real.csv");
  //  itiK.GoToBegin();
  //  for(; !itiK.IsAtEnd(); ++itiK)
  //  {
  //    out9 << itiK.Get() << std::endl;
  //  }
  //  out9.close();
  //  // DEBUG

  // square kernel and multiply with finterp
  itiK.GoToBegin();
  for (; !itiK.IsAtEnd(); ++itiK)
  {
    itiK.Set(itiK.Get() * itiK.Get() * finterp);
  }

  //  /////////////////////////////////////////////////////////////////////////////////////////////
  //  // DEBUG
  //  std::ofstream out2("kernel_sq_real.csv");
  //  itiK.GoToBegin();
  //  for(; !itiK.IsAtEnd(); ++itiK)
  //  {
  //    out2 << itiK.Get() << std::endl;
  //  }
  //  out2.close();
  //  // DEBUG

  // FFT kernel
  typename FFTType::Pointer fftK2 = FFTType::New();
  fftK2->SetInput(ifftK->GetOutput());
  fftK2->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
  fftK2->Update();
  this->m_KernelFFT = fftK2->GetOutput();

  //  /////////////////////////////////////////////////////////////////////////////////////////////
  //  // DEBUG
  //  std::ofstream out3("kernel_sq_fourier.csv");
  //  IteratorType itK2(this->m_KernelFFT, this->m_KernelFFT->GetLargestPossibleRegion() );
  //  itK2.GoToBegin();
  //  for(; !itK2.IsAtEnd(); ++itK2)
  //  {
  //    out3 << itK2.Get().real() << " " << itK2.Get().imag() << std::endl;
  //  }
  //  out3.close();
  //  // DEBUG

  // Replicate and window if required
  if (this->GetHannCutFrequencyY() > 0.)
  {
    std::cerr << "HannY not implemented for variance image filter." << std::endl;

    size.Fill(1);
    size[0] = this->m_KernelFFT->GetLargestPossibleRegion().GetSize(0);
    size[1] = height;

    const unsigned int ncut = itk::Math::Round<double>((height / 2 + 1) * std::min(1.0, this->GetHannCutFrequencyY()));

    this->m_KernelFFT = FFTOutputImageType::New();
    this->m_KernelFFT->SetRegions(size);
    this->m_KernelFFT->Allocate();
    this->m_KernelFFT->FillBuffer(0.);

    IteratorType itTwoDK(this->m_KernelFFT, this->m_KernelFFT->GetLargestPossibleRegion());
    for (unsigned int j = 0; j < ncut; j++)
    {
      itK.GoToBegin();
      const TFFTPrecision win(0.5 * (1 + std::cos(itk::Math::pi * j / ncut)));
      for (unsigned int i = 0; i < size[0]; ++itK, ++itTwoDK, i++)
      {
        itTwoDK.Set(win * itK.Get());
      }
    }
    itTwoDK.GoToReverseBegin();
    for (unsigned int j = 1; j < ncut; j++)
    {
      itK.GoToReverseBegin();
      const TFFTPrecision win(0.5 * (1 + std::cos(itk::Math::pi * j / ncut)));
      for (unsigned int i = 0; i < size[0]; --itK, --itTwoDK, i++)
      {
        itTwoDK.Set(win * itK.Get());
      }
    }
  }

  this->m_KernelFFT->DisconnectPipeline();

  //  //////////////////////////////////////////////////////////////////////////////////////////////
  //  // DEBUG
  //  std::ofstream out4("finterp.csv");
  //  out4 << finterp << std::endl;
  //  out4 << ratiogcg2 << std::endl;
  //  out4.close();
  //  std::cerr << std::endl << "finterp = " << finterp << std::endl;
  //  std::cerr << "ratiogcg2 = " << ratiogcg2 << std::endl;
  //  std::cerr << "Aborting here since filters were alrady written out." << std::endl;
  //  exit(9);
  //  // DEBUG
}

} // namespace pct
#endif
