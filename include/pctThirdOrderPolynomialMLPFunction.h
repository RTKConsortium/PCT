#ifndef __pctThirdOrderPolynomialMLPFunction_h
#define __pctThirdOrderPolynomialMLPFunction_h

#include "pctMostLikelyPathFunction.h"

namespace pct
{

/** \class ThirdOrderPolynomialMLPFunction
 * \brief Fit a third order polynomial given input and output direction and position.
 *
 * \ingroup Functions PCT
 */
template <class TCoordRep = double>
class ITK_EXPORT ThirdOrderPolynomialMLPFunction : public MostLikelyPathFunction<TCoordRep>
{
public:
  /** Standard class typedefs. */
  using Self = ThirdOrderPolynomialMLPFunction;
  using Superclass = MostLikelyPathFunction<TCoordRep>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Useful defines. */
  using VectorType = typename Superclass::VectorType;

  /** Init the mlp parameters from the input and output directions and positions. */
  virtual void
  Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut) override;

  /** Evaluate the coordinates (x,y) at depth z. */
  virtual void
  Evaluate(const TCoordRep z, TCoordRep & x, TCoordRep & y, TCoordRep & dx, TCoordRep & dy) override;

protected:
  /// Constructor
  ThirdOrderPolynomialMLPFunction() {}

  /// Destructor
  ~ThirdOrderPolynomialMLPFunction() {}

private:
  ThirdOrderPolynomialMLPFunction(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  // Polynomial parameters in each direction
  TCoordRep ax, bx, cx, dx;
  TCoordRep ay, by, cy, dy;
  TCoordRep zoffset;
};

} // namespace pct

#include "pctThirdOrderPolynomialMLPFunction.hxx"

#endif
