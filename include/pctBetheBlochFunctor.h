#ifndef __pctBetheBlochFunctor_h
#define __pctBetheBlochFunctor_h

#include "pctPhysicalConstants.h"
#include <itkImage.h>

namespace pct
{

/** \class BetheBlochProtonStoppingPower
 * \ingroup PCT
 * \brief Proton stopping power according to the Bethe-Bloch equation
 * The function corresponds to:
 *  - F in [Li, MIC-NSS, 2003].
 *  - S in [Schulte, Med Phys, 2005].
 *  - S in [Penfold, Med Phys, 2009].
 * Note that they have all an error. r must be squared in the equation.
 * The CLHEP system of units is used.
 */

namespace Functor
{

template <class TInput, class TOutput>
class BetheBlochProtonStoppingPower
{
public:
  BetheBlochProtonStoppingPower() = default;
  ~BetheBlochProtonStoppingPower() {}

  TOutput
  GetValue(const TInput e, const double I) const
  {
    /** Physical constants */
    static const double K = 4. * CLHEP::pi * CLHEP::classic_electr_radius * CLHEP::classic_electr_radius *
                            CLHEP::electron_mass_c2 * 3.343e+23 / CLHEP::cm3;

    TOutput betasq = CLHEP::proton_mass_c2 / (e + CLHEP::proton_mass_c2);
    betasq = 1. - betasq * betasq;
    return K * (log(2. * CLHEP::electron_mass_c2 / I * betasq / (1 - betasq)) - betasq) / betasq;
  }

  TOutput
  GetLowEnergyLimit()
  {
    return 1 * CLHEP::keV; // Arbitrary
  }
};

/** \class IntegratedBetheBlochProtonStoppingPowerInverse
 * \ingroup PCT
 * \brief Numerical for integral used in proton CT.
 */
template <class TInput, class TOutput>
class IntegratedBetheBlochProtonStoppingPowerInverse
{
public:
  IntegratedBetheBlochProtonStoppingPowerInverse(const double I,
                                                 const double maxEnergy,
                                                 const double binSize = 1. * CLHEP::keV)
    : m_BinSize(binSize)
  {
    unsigned int lowBinLimit, numberOfBins;
    lowBinLimit = itk::Math::Ceil<unsigned int, double>(m_S.GetLowEnergyLimit() / m_BinSize);
    numberOfBins = itk::Math::Ceil<unsigned int, double>(maxEnergy / m_BinSize);
    m_LUT.resize(numberOfBins);
    // Create lookup table for integer values of energy
    for (unsigned int i = 0; i < lowBinLimit; i++)
    {
      m_LUT[i] = 0.;
    }
    for (unsigned int i = lowBinLimit; i < numberOfBins; i++)
    {
      m_LUT[i] = m_LUT[i - 1] + binSize / m_S.GetValue(TOutput(i) * binSize, I);
      // Create inverse lut, i.e., get energy from length in water
      for (unsigned j = m_Length.size(); j < unsigned(m_LUT[i] * CLHEP::mm) + 1; j++)
      {
        m_Length.push_back(binSize * (i + double(j - m_LUT[i] * CLHEP::mm) / ((m_LUT[i] - m_LUT[i - 1]) * CLHEP::mm)));
      }
    }
  }
  ~IntegratedBetheBlochProtonStoppingPowerInverse() {}

  /** Get the integral from 0. to e. */
  TOutput
  GetValue(const TInput e) const
  {
    // SR: linear interpolation instead?
    return m_LUT[itk::Math::Round<int, TInput>(e / m_BinSize)];
  }

  /** Get the integral from e1 to e2. */
  TOutput
  GetValue(const TInput e1, const TInput e2) const
  {
    // SR: linear interpolation instead?
    return m_LUT[itk::Math::Round<int, TInput>(e2 / m_BinSize)] - m_LUT[itk::Math::Round<int, TInput>(e1 / m_BinSize)];
  }

  /** Get the energy required to traverse length l in water. */
  TOutput
  GetEnergy(const TInput l) const
  {
    // SR: linear interpolation instead?
    return m_Length[itk::Math::Round<int, TInput>(l / CLHEP::mm)];
  }

  /** Get the residual energy after length l in water with initial energy e0. */
  TOutput
  GetEnergy(const TInput l, const TInput e0) const
  {
    // SR: linear interpolation instead?
    return GetEnergy(GetValue(e0) - l);
  }

private:
  BetheBlochProtonStoppingPower<TOutput, TOutput> m_S;
  std::vector<TOutput>                            m_Length;

  double               m_BinSize;
  std::vector<TOutput> m_LUT;
};
} // end namespace Functor
} // end namespace pct

#endif
