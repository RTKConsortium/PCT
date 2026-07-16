/*=========================================================================
 *
 *  Copyright PCT Consortium
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

#ifndef pctGgoFunctions_h
#define pctGgoFunctions_h

#include "rtkMacro.h"
#include "rtkConstantImageSource.h"
#include "rtkIOFactories.h"
#include "rtkProjectionsReader.h"
#include <itkRegularExpressionSeriesFileNames.h>
#include <itksys/RegularExpression.hxx>

namespace pct
{

/** \brief Read a stack of 2D projections from gengetopt specifications.
 *
 * This function returns the file names of a projection series from command
 * line options stored in ggo struct.
 * The required options in the ggo struct are:
 *     - verbose
 *     - path: path containing projections
 *     - regexp: regular expression to select projection files in path
 *     - nsort: boolean to (des-)activate the numeric sort for expression matches
 *     - submatch: index of the submatch that will be used to sort matches
 *
 * \author Simon Rit
 *
 * \ingroup PCT
 */
template <class TArgsInfo>
std::vector<std::string>
GetProjectionsFileNamesFromGgo(const TArgsInfo & args_info)
{
  auto names = itk::RegularExpressionSeriesFileNames::New();
  names->SetDirectory(args_info.path_arg);
  names->SetNumericSort(args_info.nsort_flag);
  names->SetRegularExpression(args_info.regexp_arg);
  names->SetSubMatch(args_info.submatch_arg);

  if (args_info.verbose_flag)
    std::cout << "Regular expression matches " << names->GetFileNames().size() << " file(s)..." << std::endl;

  if (args_info.submatch_given)
  {
    itksys::RegularExpression reg;
    if (!reg.compile(args_info.regexp_arg))
    {
      itkGenericExceptionMacro(<< "Error compiling regular expression " << args_info.regexp_arg);
    }

    for (const std::string & name : names->GetFileNames())
    {
      reg.find(name);
      if (reg.match(args_info.submatch_arg) == std::string(""))
      {
        itkGenericExceptionMacro(<< "Cannot find submatch " << args_info.submatch_arg << " in " << name
                                 << " from regular expression " << args_info.regexp_arg);
      }
    }
  }

  std::vector<std::string> fileNames = names->GetFileNames();
  rtk::RegisterIOFactories();
  std::vector<size_t> idxtopop;
  size_t              i = 0;
  for (const auto & fn : fileNames)
  {
    itk::ImageIOBase::Pointer imageio =
      itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::IOFileModeEnum::ReadMode);

    if (imageio.IsNull())
    {
      std::cerr << "Ignoring file: " << fn << "\n";
      idxtopop.push_back(i);
    }
    i++;
  }
  std::reverse(idxtopop.begin(), idxtopop.end());
  for (const auto & id : idxtopop)
  {
    fileNames.erase(fileNames.begin() + id);
  }

  return fileNames;
}

template <class TProjectionsReaderType, class TArgsInfo>
void
SetProjectionsReaderFromGgo(TProjectionsReaderType * reader, const TArgsInfo & args_info)
{
  const std::vector<std::string> fileNames = GetProjectionsFileNamesFromGgo(args_info);

  const unsigned int Dimension = TProjectionsReaderType::OutputImageType::GetImageDimension();

  typename TProjectionsReaderType::OutputImageDirectionType direction;
  if (args_info.newdirection_given)
  {
    direction.Fill(args_info.newdirection_arg[0]);
    for (unsigned int i = 0; i < args_info.newdirection_given; i++)
      direction[i / Dimension][i % Dimension] = args_info.newdirection_arg[i];
    reader->SetDirection(direction);
  }

  typename TProjectionsReaderType::OutputImageSpacingType spacing;
  if (args_info.newspacing_given)
  {
    spacing.Fill(args_info.newspacing_arg[0]);
    for (unsigned int i = 0; i < args_info.newspacing_given; i++)
      spacing[i] = args_info.newspacing_arg[i];
    reader->SetSpacing(spacing);
  }

  typename TProjectionsReaderType::OutputImagePointType origin;
  if (args_info.neworigin_given)
  {
    origin.Fill(args_info.neworigin_arg[0]);
    for (unsigned int i = 0; i < args_info.neworigin_given; i++)
      origin[i] = args_info.neworigin_arg[i];
    reader->SetOrigin(origin);
  }

  auto lowerCrop = itk::MakeFilled<typename TProjectionsReaderType::OutputImageSizeType>(0);
  for (unsigned int i = 0; i < args_info.lowercrop_given; i++)
    lowerCrop[i] = args_info.lowercrop_arg[i];
  if (args_info.lowercrop_given)
    reader->SetLowerBoundaryCropSize(lowerCrop);

  auto upperCrop = itk::MakeFilled<typename TProjectionsReaderType::OutputImageSizeType>(0);
  for (unsigned int i = 0; i < args_info.uppercrop_given; i++)
    upperCrop[i] = args_info.uppercrop_arg[i];
  if (args_info.uppercrop_given)
    reader->SetUpperBoundaryCropSize(upperCrop);

  typename TProjectionsReaderType::MedianRadiusType medianRadius{};
  for (unsigned int i = 0; i < args_info.radius_given; i++)
    medianRadius[i] = args_info.radius_arg[i];
  reader->SetMedianRadius(medianRadius);
  if (args_info.multiplier_given)
    reader->SetConditionalMedianThresholdMultiplier(args_info.multiplier_arg);

  typename TProjectionsReaderType::ShrinkFactorsType binFactors;
  binFactors.Fill(1);
  for (unsigned int i = 0; i < args_info.binning_given; i++)
    binFactors[i] = args_info.binning_arg[i];
  reader->SetShrinkFactors(binFactors);

  if (args_info.wpc_given)
  {
    std::vector<double> coeffs;
    coeffs.assign(args_info.wpc_arg, args_info.wpc_arg + args_info.wpc_given);
    reader->SetWaterPrecorrectionCoefficients(coeffs);
  }

  reader->SetFileNames(fileNames);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(reader->UpdateOutputInformation());
}

} // namespace pct

#endif // pctGgoFunctions_h
