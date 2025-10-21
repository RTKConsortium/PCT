/*=========================================================================
 *
 *  Copyright Centre National de la Recherche Scientifique
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "pctPolynomialMLPFunction.h"

#include "itkTestingMacros.h"

#include <vector>
#include <iostream>

int
pctPolynomialMLPFunctionTest(int, char **)
{
  auto mlp = pct::PolynomialMLPFunction::New();

  // Simple geometry for Init()
  using VectorType = pct::MostLikelyPathFunction<double>::VectorType;
  VectorType posIn = itk::MakeVector(0.0, 0.0, 0.0);
  VectorType posOut = itk::MakeVector(0.0, 0.0, 100.0);
  VectorType dirIn = itk::MakeVector(0.0, 0.0, 1.0);
  VectorType dirOut = itk::MakeVector(0.0, 0.0, 1.0);

  // Populate built-in coefficients
  mlp->Init(posIn, posOut, dirIn, dirOut);

  // Default degree must be 5
  std::vector<double> coeffs = mlp->GetCoefficients();
  if (coeffs.size() != 6)
  {
    std::cerr << "Default coefficients size != 6 (degree 5). Got " << coeffs.size() << std::endl;
    return EXIT_FAILURE;
  }

  // Set degree to 2 and re-Init
  mlp->SetPolynomialDegree(2);
  mlp->Init(posIn, posOut, dirIn, dirOut);

  std::vector<double> u = { -50.0, 0.0, 50.0 };
  std::vector<double> x, y;

  // Test MostLikelyPathFunction::Evaluate
  mlp->Evaluate(u, x, y);

  // Basic sanity checks on Evaluate output: correct sizes
  if (x.size() != u.size() || y.size() != u.size())
  {
    std::cerr << "Evaluate returned vectors of wrong size." << std::endl;
    return EXIT_FAILURE;
  }

  // Custom coefficients (degree inferred from vector length). Use degree = 3 (4 coefficients).
  std::vector<double> custom = { 1.0e-06, 2.0e-07, 3.0e-08, 4.0e-09 };
  mlp->SetCoefficients(custom);
  mlp->Init(posIn, posOut, dirIn, dirOut);
  // Verify coefficients round-trip through GetCoefficients
  std::vector<double> got = mlp->GetCoefficients();
  if (got.size() != custom.size())
  {
    std::cerr << "GetCoefficients size mismatch: " << got.size() << " != " << custom.size() << std::endl;
    return EXIT_FAILURE;
  }
  for (unsigned int i = 0; i < got.size(); ++i)
  {
    if (got[i] != custom[i])
    {
      std::cerr << "Coefficient mismatch at " << i << ": " << got[i] << " != " << custom[i] << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "pctPolynomialMLPFunctionTest PASSED" << std::endl;
  return EXIT_SUCCESS;
}
