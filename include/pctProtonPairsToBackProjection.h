#ifndef __pctProtonPairsToBackProjection_h
#define __pctProtonPairsToBackProjection_h

#include "rtkConfiguration.h"
#include "pctBetheBlochFunctor.h"

#include <rtkQuadricShape.h>
#include <rtkThreeDCircularProjectionGeometry.h>
#include <itkInPlaceImageFilter.h>
#include <mutex>

namespace pct
{

template <class TInputImage, class TOutputImage>
class ITK_TEMPLATE_EXPORT ProtonPairsToBackProjection : public itk::InPlaceImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  using Self = ProtonPairsToBackProjection;
  using Superclass = itk::InPlaceImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using FileNamesContainer = std::vector<std::string>;

  using ProtonPairsPixelType = itk::Vector<float, 3>;
  using ProtonPairsImageType = itk::Image<ProtonPairsPixelType, 2>;
  using ProtonPairsImagePointer = ProtonPairsImageType::Pointer;

  using CountImageType = typename itk::Image<unsigned int, TInputImage::ImageDimension>;
  using CountImagePointer = typename CountImageType::Pointer;

  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputImageRegionType = typename OutputImageType::RegionType;

  using RQIType = rtk::QuadricShape;

  using GeometryType = rtk::ThreeDCircularProjectionGeometry;
  using GeometryPointer = typename GeometryType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(ProtonPairsToBackProjection);

  /** Set the vector of strings that contains the file names. Files
   * are processed in sequential order. */
  void
  SetProtonPairsFileNames(const FileNamesContainer & name)
  {
    if (m_ProtonPairsFileNames != name)
    {
      m_ProtonPairsFileNames = name;
      this->Modified();
    }
  }
  const FileNamesContainer &
  GetProtonPairsFileNames() const
  {
    return m_ProtonPairsFileNames;
  }

  /** Get/Set the most likely path type. Can be "schulte" or "polynomial" */
  itkGetMacro(MostLikelyPathType, std::string);
  itkSetMacro(MostLikelyPathType, std::string);

  itkGetMacro(MostLikelyPathPolynomialDegree, int);
  itkSetMacro(MostLikelyPathPolynomialDegree, int);

  /** Get/Set the boundaries of the object. */
  itkGetMacro(QuadricIn, RQIType::Pointer);
  itkSetMacro(QuadricIn, RQIType::Pointer);
  itkGetMacro(QuadricOut, RQIType::Pointer);
  itkSetMacro(QuadricOut, RQIType::Pointer);

  /** Get/Set the count of proton pairs per pixel. */
  itkGetMacro(Counts, CountImagePointer);
  itkSetMacro(Counts, CountImagePointer);

  /** Get/Set the ionization potential used in the Bethe-Bloch equation. */
  itkGetMacro(IonizationPotential, double);
  itkSetMacro(IonizationPotential, double);

  /** Get / Set the object pointer to projection geometry */
  itkGetMacro(Geometry, GeometryPointer);
  itkSetMacro(Geometry, GeometryPointer);

  /** Get / Set the boolean desactivating rotation to bin in coordinate orientation. Default is off. */
  itkGetMacro(DisableRotation, bool);
  itkSetMacro(DisableRotation, bool);
  itkBooleanMacro(DisableRotation);

protected:
  ProtonPairsToBackProjection();
  virtual ~ProtonPairsToBackProjection() {}

  virtual void
  BeforeThreadedGenerateData() override;
  virtual void
  GenerateData() override;
  virtual void
  AfterThreadedGenerateData() override;

  /** The two inputs should not be in the same space so there is nothing
   * to verify. */
  virtual void
  VerifyInputInformation() const override
  {}

private:
  ProtonPairsToBackProjection(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  /** A list of filenames to be processed. */
  FileNamesContainer m_ProtonPairsFileNames;

  std::string m_MostLikelyPathType;
  int         m_MostLikelyPathPolynomialDegree;

  /** Count event in each thread */
  CountImagePointer m_Counts;

  /** The two quadric functions defining the object support. */
  RQIType::Pointer m_QuadricIn;
  RQIType::Pointer m_QuadricOut;

  /** Ionization potential used in the Bethe Bloch equation */
  double m_IonizationPotential;

  /** The functor to convert energy loss to attenuation */
  Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double> * m_ConvFunc;

  ProtonPairsImageType::Pointer m_ProtonPairs;

  /** RTK geometry object */
  GeometryPointer m_Geometry;

  /** Disable rotation to bin in coordinate orientation. Default is off. */
  bool m_DisableRotation = false;

  std::mutex m_Mutex;
};

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#  include "pctProtonPairsToBackProjection.hxx"
#endif

#endif
