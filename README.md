# PCT

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/RTKConsortium/PCT/blob/master/LICENSE)
[![][gha-img]][gha-link]

PCT (Proton Computed Tomography) is a toolkit used to process proton CT data and reconstruct proton stopping power maps. PCT is written both in C++ and Python, and is designed to be used either as a code library, or with command-line applications.

## Links

* [Documentation](https://proton-ct.readthedocs.io)
* [Issue tracking](https://github.com/RTKConsortium/PCT/issues)
* [ITK](https;//itk.org)

## Usage

Usage of each PCT application can be described using the `--help`/`-h` option. For instance, running
```bash
pctfdk --help
```
displays the help for `pctfdk`.

Reconstruction typically involves the following steps:
- `pctpairprotons` in order to arrange ROOT data in a format described [here](pct_format.md).
- `pctpaircuts` in order to remove nuclear collisions.
- `pctbinning` in order to compute the distance-driven binning as described [here]( https://doi.org/10.1118/1.4789589).
- `pctfdk` in order to reconstruct the data generated in the previous step using distance-driven FDK (as described [here](https://doi.org/10.1118/1.4789589)).

## Copyright Centre National de la Recherche Scientifique

Licensed under the Apache License, Version 2.0 (the "[License](https://www.apache.org/licenses/LICENSE-2.0.txt)"); you may not use this file except in compliance with the [License](https://www.apache.org/licenses/LICENSE-2.0.txt).

Unless required by applicable law or agreed to in writing, software distributed under the [License](https://www.apache.org/licenses/LICENSE-2.0.txt) is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the [License](https://www.apache.org/licenses/LICENSE-2.0.txt) for the specific language governing permissions and limitations under the [License](https://www.apache.org/licenses/LICENSE-2.0.txt).
