# Installation

## User version

The easiest way to install PCT and its dependencies is using the automatically generated Python wheels that can be downloaded directly on [PCT's GitHub Actions page](https://github.com/RTKConsortium/PCT/actions). Log in to GitHub, navigate to the artifact section of an action that recently finished successfully, and download the wheel that corresponds to your operating system and architecture. For example, a Linux user on an x64 machine should download the file named `LinuxWheel311_2_28-x64`. Then unzip the downloaded folder. It should contain a file with a name like `itk_pct-0.1.0-cp311-abi3-manylinux_2_28_x86_64.whl`.

The `311` in the name above correspond to the Python version, in that case Python 3.11. Python in backward-compatible, so Python 3.12 or above also be able to run the wheel. If the Python version of your computer matches this requirement, then you can install it. However, it is much advised to intall PCT in a virtual environment, created for instance unsing [`conda`](https://conda.io) or [`virtualenv`](https://virtualenv.pypa.io). Older versions of Python are compiled on the master branch only.

Once a proper Python version is set up, PCT can be installed using
```bash
pip install <path to the wheel file>
```

You can check that PCT was properly installed by opening a Python interpreter, then trying to import PCT using
```python
from itk import PCT
```
If no error appears, then PCT was successfully installed.

## Developer version

If you need to edit the source code of PCT, or tune the compilation process, then you need to manually compile PCT. PCT depends on [ITK](https://itk.org) and [RTK](https://openrtk.org), the first step is therefore to download and compile them.

### Compiling ITK and RTK

RTK comes bundled as an ITK module, so it is enough to download ITK and compile it with an option that indicates that RTK is also required. First, clone ITK from GitHub:
```bash
git clone git@github.com:InsightSoftwareConsortium/ITK.git
```

Then, create a compilation folder for ITK:
```
mkdir itk-build
cd itk-build
```

ITK uses [CMake](https://cmake.org) as its building tool. Below is an example of how to configure ITK, but feel free to edit the options as needed.
```bash
cmake \
    -DITK_WRAP_PYTHON=ON \
    -DModule_RTK=ON \
    ../ITK
```
`ITK_WRAP_PYTHON` mean that ITK will generate and compile Python wrapping for the C++ code, which is also supported by PCT. This option massively increases compilation time, so if Python wrappings are not needed, this option can be disabled. `Module_RTK` enables the compilation of RTK, which is required for PCT.

Once configuration finishes, the compilation is started using
```
cmake --build .
```
Depending on your machine and whether `ITK_WRAP_PYTHON` is on, compiling ITK can take from few minutes up to several hours.

### Compiling PCT

The source code of PCT can be downloaded from it [git repository](https://github.com/RTKConsortium/PCT) hosted at GitHub. The command to clone the repository is:
```bash
git clone https://github.com/RTKConsortium/PCT.git
```

Then, create a build directory for PCT:
```bash
mkdir pct-build
cd pct-build
```

Similar to ITK, configuration and compilation is handled via [CMake](https://cmake.org/). Configuration can be achieved using
```bash
cmake \
    -DITK_DIR=<path to ITK build folder> \
    -DPCT_BUILD_APPLICATIONS=ON \
    ../PCT
```
Once configuration completes, the compilation can be achieved using
```bash
cmake --build .
```

Optionally, PCT can be added to the user's `$PATH` variable using for instance
```bash
export PATH=<path to ITK build folder>/Wrapping/Generators/Python/itk/:${PATH}  # where PCT was built
```
in the user's `.bashrc` file. This allows to run PCT applications from anywhere on the computer.

## Optional dependencies

### GATE

[GATE](http://www.opengatecollaboration.org/) is an open source software for Monte Carlo simulation of emission tomography, computed tomography, optical imaging and radiotherapy experiments, based on [Geant4](https://cern.ch/geant4). Some features of PCT require to run GATE. [Instructions on how to install GATE can be found in the GATE documentation](https://opengate-python.readthedocs.io/en/master/user_guide/user_guide_installation.html), but the simplest way is to install the `opengate` package:
```bash
pip install opengate
```
An example of proton CT GATE simulation can be found in [`gate/protonct.py`](https://github.com/RTKConsortium/PCT/blob/master/gate/protonct.py).
