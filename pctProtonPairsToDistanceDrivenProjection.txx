#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>

#include "pctThirdOrderPolynomialMLPFunction.h"
#include "pctSchulteMLPFunction.h"
#include "pctEnergyStragglingFunctor.h"

namespace pct
{

template <class TInputImage, class TOutputImage>
ProtonPairsToDistanceDrivenProjection<TInputImage, TOutputImage>
::ProtonPairsToDistanceDrivenProjection():m_Robust(false),m_ComputeScattering(false)
{
  this->DynamicMultiThreadingOff();
  this->SetNumberOfWorkUnits( itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads() );
}

template <class TInputImage, class TOutputImage>
void
ProtonPairsToDistanceDrivenProjection<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  m_Outputs.resize( this->GetNumberOfWorkUnits() );
  m_Counts.resize( this->GetNumberOfWorkUnits() );
  if(m_ComputeScattering)
    {
    m_Angles.resize( this->GetNumberOfWorkUnits() );
    m_AnglesVectors.resize( this->GetInput()->GetLargestPossibleRegion().GetNumberOfPixels() );
    m_AnglesSq.resize( this->GetNumberOfWorkUnits() );
    }
  if(m_SigmaMap)
    {
    m_Sigmas.resize( this->GetNumberOfWorkUnits() );
    }

  if(m_QuadricOut.GetPointer()==NULL)
    m_QuadricOut = m_QuadricIn;
  m_ConvFunc = new Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double>(m_IonizationPotential, 600.*CLHEP::MeV, 0.1*CLHEP::keV);

  // Read pairs
  typedef itk::ImageFileReader< ProtonPairsImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( m_ProtonPairsFileName );
  reader->Update();
  m_ProtonPairs = reader->GetOutput();
}

template <class TInputImage, class TOutputImage>
void
ProtonPairsToDistanceDrivenProjection<TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& itkNotUsed(outputRegionForThread), rtk::ThreadIdType threadId)
{
  // Create MLP depending on type
  pct::SchulteMLPFunction::Pointer mlp = pct::SchulteMLPFunction::New();

  // Create thread image and corresponding stack to count events
  m_Counts[threadId] = CountImageType::New();
  m_Counts[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
  m_Counts[threadId]->Allocate();
  m_Counts[threadId]->FillBuffer(0);

  if( m_ComputeScattering && (!m_Robust || threadId==0) )
    {
    m_Angles[threadId] = AngleImageType::New();
    m_Angles[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_Angles[threadId]->Allocate();
    m_Angles[threadId]->FillBuffer(0);

    m_AnglesSq[threadId] = AngleImageType::New();
    m_AnglesSq[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_AnglesSq[threadId]->Allocate();
    m_AnglesSq[threadId]->FillBuffer(0);
    }

  if( m_SigmaMap )
    {
    m_Sigmas[threadId] = OutputImageType::New();
    m_Sigmas[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_Sigmas[threadId]->Allocate();
    m_Sigmas[threadId]->FillBuffer(0);
    }

  if(threadId==0)
    {
    m_Outputs[0] = this->GetOutput();
    m_Count = m_Counts[0];
    if(m_ComputeScattering)
      {
      m_Angle = m_Angles[0];
      m_AngleSq = m_AnglesSq[0];
      }
    if(m_SigmaMap)
      {
      m_Sigma = m_Sigmas[0];
      }
    }
  else
    {
    m_Outputs[threadId] = OutputImageType::New();
    m_Outputs[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_Outputs[threadId]->Allocate();
    }
  m_Outputs[threadId]->FillBuffer(0.);

  size_t nprotons = m_ProtonPairs->GetLargestPossibleRegion().GetSize()[1];
  ProtonPairsImageType::RegionType region = m_ProtonPairs->GetLargestPossibleRegion();
  region.SetIndex(1, threadId*nprotons/this->GetNumberOfWorkUnits());
  region.SetSize(1, std::min((unsigned long)nprotons/this->GetNumberOfWorkUnits(), nprotons-region.GetIndex(1)));

  // Image information constants
  const typename OutputImageType::SizeType    imgSize    = this->GetInput()->GetBufferedRegion().GetSize();
  const typename OutputImageType::PointType   imgOrigin  = this->GetInput()->GetOrigin();
  const typename OutputImageType::SpacingType imgSpacing = this->GetInput()->GetSpacing();
  const unsigned long npixelsPerSlice = imgSize[0] * imgSize[1];

  typename OutputImageType::PixelType *imgData = m_Outputs[threadId]->GetBufferPointer();
  float *imgCountData = m_Counts[threadId]->GetBufferPointer();
  float *imgAngleData = NULL, *imgAngleSqData = NULL;
  float *imgSigmaData = NULL;
  if(m_ComputeScattering && !m_Robust)
    {
    imgAngleData = m_Angles[threadId]->GetBufferPointer();
    imgAngleSqData = m_AnglesSq[threadId]->GetBufferPointer();
    }
  if(m_SigmaMap)
    {
    imgSigmaData = m_Sigmas[threadId]->GetBufferPointer();
    }

  itk::Vector<float, 3> imgSpacingInv;
  for(unsigned int i=0; i<3; i++)
    imgSpacingInv[i] = 1./imgSpacing[i];

  // Corrections
  typedef itk::Vector<double,3> VectorType;

  // Create zmm and magnitude lut (look up table)
  itk::ImageRegionIterator<ProtonPairsImageType> it(m_ProtonPairs, region);
  std::vector<double> zmm(imgSize[2]);
  std::vector<double> zmag(imgSize[2]);
  ++it;
  const double zPlaneOutInMM = it.Get()[2];
  --it;
  for(unsigned int i=0; i<imgSize[2]; i++)
    {
    zmm[i] = i*imgSpacing[2]+imgOrigin[2];
    zmag[i] = (m_SourceDistance==0.)?1:(zPlaneOutInMM-m_SourceDistance)/(zmm[i]-m_SourceDistance);
    }

  GeometryType::ThreeDHomogeneousMatrixType rotMat;
  rotMat = m_Geometry->GetRotationMatrices()[m_Index].GetInverse();

// Process pairs
  while(!it.IsAtEnd() ) //
    {
    if(threadId==0 && it.GetIndex()[1]%10000==0)
      {
      std::cout << '\r'
                << it.GetIndex()[1] << " pairs of protons processed ("
                << 100*it.GetIndex()[1]/region.GetSize(1) << "%) in thread 1"
                << std::flush;
      }
    VectorType pIn = it.Get();
    ++it;
    VectorType pOut = it.Get();
    ++it;
    VectorType dIn = it.Get();
    ++it;
    VectorType dOut = it.Get();
    ++it;

    double anglex = 0., angley = 0.;
    if( m_ComputeScattering )
      {
      typedef itk::Vector<double, 2> VectorTwoDType;

      VectorTwoDType dInX, dInY, dOutX, dOutY;
      dInX[0] = dIn[0];
      dInX[1] = dIn[2];
      dInY[0] = dIn[1];
      dInY[1] = dIn[2];
      dOutX[0] = dOut[0];
      dOutX[1] = dOut[2];
      dOutY[0] = dOut[1];
      dOutY[1] = dOut[2];

      angley = std::acos( std::min(1.,dInY*dOutY / ( dInY.GetNorm() * dOutY.GetNorm() ) ) );
      anglex = std::acos( std::min(1.,dInX*dOutX / ( dInX.GetNorm() * dOutX.GetNorm() ) ) );
      }

    if(pIn[2] > pOut[2])
      {
      itkGenericExceptionMacro("Required condition pIn[2] > pOut[2] is not met, check coordinate system.");
      }
    if(dIn[2] < 0.)
      {
      itkGenericExceptionMacro("The code assumes that protons move in positive z.");
      }

    const double eIn = it.Get()[0];
    const double eOut = it.Get()[1];

    double value = 0.;
    if(eIn==0.)
      value = eOut; // Directly read WEPL
    else
      {
      value = m_ConvFunc->GetValue(eOut, eIn); // convert to WEPL
      }
    ++it;

    VectorType nucInfo(0.);
    if(it.GetIndex()[0] != 0)
      {
      nucInfo = it.Get();
      ++it;
      }
    std::vector<typename OutputImageType::OffsetValueType> offsets;
    
    VectorType pRotIn (0.);
    VectorType dRotIn (0.);
    VectorType pRotOut(0.);
    VectorType dRotOut (0.);
    for(unsigned int i=0; i<3; i++)
      {
      for(unsigned int j=0; j<3; j++)
        {
        pRotIn[i] += rotMat[i][j] * pIn[j];
        dRotIn[i] += rotMat[i][j] * dIn[j];
        pRotOut[i] += rotMat[i][j] * pOut[j];
        dRotOut[i] += rotMat[i][j] * dOut[j];
        }
      }

    // Move straight to entrance and exit shapes
    VectorType pSIn  = pIn;
    VectorType pSOut = pOut;
    double nearDistIn, nearDistOut, farDistIn, farDistOut;
    double dHullIn, dHullOut;
    bool SLP=0;
    if(m_QuadricIn.GetPointer()!=NULL)
      {
      if(m_QuadricIn->IsIntersectedByRay(pRotIn,dRotIn,nearDistIn,farDistIn) &&
         m_QuadricOut->IsIntersectedByRay(pRotOut,dRotOut,nearDistOut,farDistOut))
        {
        pSIn  = pIn  + dIn  * nearDistIn;
        dHullIn=std::abs(nearDistIn);
        if(pSIn[2]<pIn[2]  || pSIn[2]>pOut[2])
        {
          pSIn  = pIn  + dIn  * farDistIn;
          dHullIn=std::abs(farDistIn);
        }
        pSOut = pOut + dOut * nearDistOut;
        dHullOut=std::abs(nearDistOut);
        if(pSOut[2]<pIn[2] || pSOut[2]>pOut[2])
        {
          pSOut = pOut + dOut * farDistOut;
          dHullOut=std::abs(farDistOut);
        }

        }
      else
        {
        SLP = 1;
        }
      }

    // Normalize direction with respect to z
    dIn[0] /= dIn[2];
    dIn[1] /= dIn[2];
    //dIn[2] = 1.; SR: implicit in the following
    dOut[0] /= dOut[2];
    dOut[1] /= dOut[2];
    //dOut[2] = 1.; SR: implicit in the following

    itk::Matrix<double, 2, 2> mlp_error;
    mlp_error.Fill(0);
    // Init MLP before mm to voxel conversion

    double xSIn = pSIn[0];
    double ySIn = pSIn[1];
    double xSOut = pSOut[0];
    double ySOut = pSOut[1];
    double dxSIn, dySIn, dxSOut, dySOut;
    VectorType dSIn = dIn;
    VectorType dSOut = dOut;
    mlp->Init(pSIn, pSOut, dIn, dOut);
    if(m_MLPKrah && !SLP)
      {
      mlp->SetTrackerInfo(std::abs(dHullIn),std::abs(dHullOut),10,0.065817931,0.005);
      mlp->EvaluateErrorWithTrackerUncertainty(pSIn[2],eIn, eOut, mlp_error,xSIn,ySIn, dxSIn, dySIn);
      mlp->EvaluateErrorWithTrackerUncertainty(pSOut[2],eIn, eOut, mlp_error,xSOut,ySOut,dxSOut,dySOut);
      dSIn[0] = std::tan(dxSIn);
      dSIn[1] = std::tan(dySIn);
      dSOut[0] = std::tan(dxSOut);
      dSOut[1] = std::tan(dySOut);
      pSIn[0] = xSIn;
      pSIn[1] = ySIn;
      pSOut[0] = xSOut;
      pSOut[1] = ySOut;
      }

    for(unsigned int k=0; k<imgSize[2]; k++)
      {
      double std_error = 0;
      double xx, yy;
      double dummy_xx, dummy_yy, dummy_dxx, dummy_dyy;
      const double dk = zmm[k];
      if(SLP)
        {
        const double z = (dk-pIn[2]);
        xx = pIn[0]+z*dIn[0];
        yy = pIn[1]+z*dIn[1];
        }
      else{
      if(dk<=pSIn[2]) //before entrance
        {
        const double z = (dk-pSIn[2]);
        xx = pSIn[0]+z*dSIn[0];
        yy = pSIn[1]+z*dSIn[1];
        if(m_MLPKrah)
          {
          mlp->EvaluateErrorWithTrackerUncertainty(zmm[k],eIn, eOut, mlp_error,dummy_xx,dummy_yy, dummy_dxx, dummy_dyy);
          std_error=mlp_error(0,0)*zmag[k]*zmag[k];
          }
        }
      else if(dk>=pSOut[2]) //after exit
        {
        const double z = (dk-pSOut[2]);
        xx = pSOut[0]+z*dSOut[0];
        yy = pSOut[1]+z*dSOut[1];
        if(m_MLPKrah)
          {
          mlp->EvaluateErrorWithTrackerUncertainty(zmm[k],eIn, eOut, mlp_error,dummy_xx,dummy_yy, dummy_dxx, dummy_dyy);
          std_error=mlp_error(0,0)*zmag[k]*zmag[k];
          }
        }
      else //MLP
        {
        if(m_MLPKrah)
          {
          mlp->EvaluateErrorWithTrackerUncertainty(zmm[k],eIn, eOut, mlp_error,xx,yy, dummy_dxx, dummy_dyy);
          std_error=mlp_error(0,0)*zmag[k]*zmag[k];
          }
        else
          {
          mlp->Evaluate(zmm[k], xx, yy);
          mlp->EvaluateError(zmm[k],mlp_error);
          std_error=mlp_error(0,0)*zmag[k]*zmag[k];
          }
        }
      }

      // Source at (0,0,args_info.source_arg), mag then to voxel
      xx = (xx*zmag[k] - imgOrigin[0]) * imgSpacingInv[0];
      yy = (yy*zmag[k] - imgOrigin[1]) * imgSpacingInv[1];

      // Lattice conversion
      const int i = itk::Math::Round<int,double>(xx);
      const int j = itk::Math::Round<int,double>(yy);
      if(i>=0 && i<(int)imgSize[0] &&
         j>=0 && j<(int)imgSize[1])
        {
        //const unsigned long idx;
        const unsigned long idx = i+j*imgSize[0]+k*npixelsPerSlice;
        if(m_WeightsCF)
          {
          offsets.push_back(i+j*imgSize[0]);
          }
        else
         {
          imgData[ idx ] += value;
          imgCountData[ idx ]++;
         }

        if(m_SigmaMap && std_error==std_error)
          {
          imgSigmaData[ idx ] += std_error;
          }

        if(m_ComputeScattering)
          {
          if(m_Robust)
            {
            m_AnglesVectorsMutex.lock();
            m_AnglesVectors[idx].push_back(anglex);
            m_AnglesVectors[idx].push_back(angley);
            m_AnglesVectorsMutex.unlock();
            }
          else
            {
            imgAngleData[ idx ] += anglex;
            imgAngleData[ idx ] += angley;
            imgAngleSqData[ idx ] += anglex*anglex;
            imgAngleSqData[ idx ] += angley*angley;
            }
          }
        }
      }

    if(m_WeightsCF)
      {
      std::vector<typename OutputImageType::OffsetValueType> channels(offsets);
      std::sort( channels.begin(), channels.end() );
      channels.erase( std::unique( channels.begin(), channels.end() ), channels.end() );
      for(unsigned int k=0; k<channels.size(); k++)
        {
        double weight =  (double) std::count(offsets.begin(),offsets.end(),channels[k])/imgSize[2];
        imgData[ channels[k] ] += value * weight * weight;
        imgCountData[ channels[k] ] += weight * weight;
        }
      }
  }

  if(threadId==0)
    {
    std::cout << '\r'
              << region.GetSize(1) << " pairs of protons processed (100%) in thread 1"
              << std::endl;
#ifdef MLP_TIMING
    mlp->PrintTiming(std::cout);
#endif
    }

}

template <class TInputImage, class TOutputImage>
void
ProtonPairsToDistanceDrivenProjection<TInputImage, TOutputImage>
::AfterThreadedGenerateData()
{
  typedef typename itk::ImageRegionIterator<TOutputImage> ImageIteratorType;
  ImageIteratorType itOut(m_Outputs[0], m_Outputs[0]->GetLargestPossibleRegion());

  typedef itk::ImageRegionIterator<CountImageType> ImageCountIteratorType;
  ImageCountIteratorType itCOut(m_Counts[0], m_Outputs[0]->GetLargestPossibleRegion());

  // Merge the projection computed in each thread to the first one
  for(unsigned int i=1; i<this->GetNumberOfWorkUnits(); i++)
    {
    if(m_Outputs[i].GetPointer() == NULL)
      continue;
    ImageIteratorType itOutThread(m_Outputs[i], m_Outputs[i]->GetLargestPossibleRegion());
    ImageCountIteratorType itCOutThread(m_Counts[i], m_Outputs[i]->GetLargestPossibleRegion());

    while(!itOut.IsAtEnd())
      {
      itOut.Set(itOut.Get()+itOutThread.Get());
      ++itOutThread;
      ++itOut;

      itCOut.Set(itCOut.Get()+itCOutThread.Get());
      ++itCOutThread;
      ++itCOut;
      }

    itOut.GoToBegin();
    itCOut.GoToBegin();
    }

  // Set count image information
  m_Count->SetSpacing( this->GetOutput()->GetSpacing() );
  m_Count->SetOrigin( this->GetOutput()->GetOrigin() );

  // Normalize eloss wepl with proton count (average)
  while(!itCOut.IsAtEnd())
    {
    if(itCOut.Get())
      itOut.Set(itOut.Get()/itCOut.Get());
    ++itOut;
    ++itCOut;
    }

  if(m_SigmaMap)
    {
    ImageIteratorType itSigmaOut(m_Sigmas[0], m_Outputs[0]->GetLargestPossibleRegion());

    // Merge the projection computed in each thread to the first one
    for(unsigned int i=1; i<this->GetNumberOfWorkUnits(); i++)
      {
      if(m_Outputs[i].GetPointer() == NULL)
        continue;
      ImageIteratorType itSigmaOutThread(m_Sigmas[i], m_Outputs[i]->GetLargestPossibleRegion());

      while(!itSigmaOut.IsAtEnd())
	{
	itSigmaOut.Set(itSigmaOut.Get()+itSigmaOutThread.Get());
	++itSigmaOutThread;
	++itSigmaOut;
	}

      itSigmaOut.GoToBegin();
      }

    // Set count image information
    m_Sigma->SetSpacing( this->GetOutput()->GetSpacing() );
    m_Sigma->SetOrigin( this->GetOutput()->GetOrigin() );

    itCOut.GoToBegin();
    // Normalize variance with proton count (average)
    while(!itCOut.IsAtEnd())
      {
      if(itCOut.Get())
      itSigmaOut.Set(itSigmaOut.Get()/itCOut.Get());
      ++itSigmaOut;
      ++itCOut;
      }
    }

  if(m_ComputeScattering)
    {
    typedef itk::ImageRegionIterator<AngleImageType> ImageAngleIteratorType;
    ImageAngleIteratorType itAngleOut(m_Angles[0], m_Outputs[0]->GetLargestPossibleRegion());

    typedef itk::ImageRegionIterator<AngleImageType> ImageAngleSqIteratorType;
    ImageAngleSqIteratorType itAngleSqOut(m_AnglesSq[0], m_Outputs[0]->GetLargestPossibleRegion());

    if(!m_Robust)
      {
      for(unsigned int i=1; i<this->GetNumberOfWorkUnits(); i++)
        {
        if(m_Outputs[i].GetPointer() == NULL)
          continue;
        ImageAngleIteratorType itAngleOutThread(m_Angles[i], m_Outputs[i]->GetLargestPossibleRegion());
        ImageAngleSqIteratorType itAngleSqOutThread(m_AnglesSq[i], m_Outputs[i]->GetLargestPossibleRegion());

        while(!itAngleOut.IsAtEnd())
          {
          itAngleOut.Set(itAngleOut.Get()+itAngleOutThread.Get());
          ++itAngleOutThread;
          ++itAngleOut;

          itAngleSqOut.Set(itAngleSqOut.Get()+itAngleSqOutThread.Get());
          ++itAngleSqOutThread;
          ++itAngleSqOut;
          }
        itAngleOut.GoToBegin();
        itAngleSqOut.GoToBegin();
        }
      }

    // Set scattering wepl image information
    m_Angle->SetSpacing( this->GetOutput()->GetSpacing() );
    m_Angle->SetOrigin( this->GetOutput()->GetOrigin() );

    itCOut.GoToBegin();
    std::vector< std::vector<float> >::iterator itAnglesVectors = m_AnglesVectors.begin();
    while(!itCOut.IsAtEnd())
      {
      if(itCOut.Get())
        {
        // Calculate angular variance (sigma2) and convert to scattering wepl
        if(m_Robust)
          {
          if(itCOut.Get()==1)
            {
            itAngleOut.Set( 0. );
            }
          else
            {
            // Angle: 38.30% (0.5 sigma) with interpolation (median is 0. and we only have positive values
            double sigmaAPos = itAnglesVectors->size()*0.3830;
            unsigned int sigmaASupPos = itk::Math::Ceil<unsigned int, double>(sigmaAPos);
            std::partial_sort(itAnglesVectors->begin(),
                              itAnglesVectors->begin()+sigmaASupPos+1,
                              itAnglesVectors->end());
            double sigmaADiff = sigmaASupPos-sigmaAPos;
            double sigma = 2.*(*(itAnglesVectors->begin()+sigmaASupPos)*(1.-sigmaADiff)+
                               *(itAnglesVectors->begin()+sigmaASupPos-1)*sigmaADiff); //x2 to get 1sigma
            itAngleOut.Set( sigma * sigma );
            }
          }
        else
          {
          double sigma2 = itAngleSqOut.Get()/itCOut.Get()/2;
          itAngleOut.Set( sigma2 );
          }
        }

      ++itCOut;
      ++itAngleOut;
      ++itAngleSqOut;
      ++itAnglesVectors;
      }
    }

  // Free images created in threads
  m_Outputs.resize( 0 );
  m_Counts.resize( 0 );
  m_Angles.resize( 0 );
  m_AnglesSq.resize( 0 );
  m_AnglesVectors.resize( 0 );
  m_Sigmas.resize( 0 );
}

}
