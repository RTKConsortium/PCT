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

#include "pctProtonPairsToDistanceDrivenProjection.h"

#include "itkTestingMacros.h"

int
pctProtonPairsToDistanceDrivenProjectionTest(int argc, char* argv[])
{
  constexpr unsigned int Dimension = 3;
  using PixelType = float;
  using ImageType = itk::Image<PixelType, Dimension>;

  using FilterType = pct::ProtonPairsToDistanceDrivenProjection<ImageType, ImageType>;
  FilterType::Pointer filter = FilterType::New();

  ITK_EXERCISE_BASIC_OBJECT_METHODS(filter, ProtonPairsToDistanceDrivenProjection, InPlaceImageFilter);

  
  return EXIT_SUCCESS;
}
