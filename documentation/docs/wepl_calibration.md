# WEPL calibration

PCT supports several calibration curve methods to convert from a given measurement (energy loss or TOF) to the corresponding WEPL.
First, the calibration curve(s) must be generated using one of the applications provided by PCT. Generating the calibration curves typically requires GATE, which can be installed following the instructions in the [installation guide](installation.md).
Then, the resulting curves can be passed to `pctpairprotons` in order to properly convert the measurements to WEPL when forming proton pairs. The resulting pairs are directly associated to the corresponding WEPL following the convention explained in the [PCT data format](pct_format.md), that is that $e_\text{in}=0$ and $e_\text{out}=\text{WEPL}$.

## Single 1D calibration curve

PCT provides the application `pctweplfit` that generates 1D calibration curves from energy loss or TOF to WEPL, following the technique described in [Coussat et al, Fully 3D, 2025](https://hal.science/hal-05375602). The resulting polynomial fit can then be used when pairing protons with `pctpairprotons` to convert the energy loss or the TOF directly to WEPL.

Below is an example of how to run `pctweplfit`:
```bash
pctweplfit \
    -o output \
    --savefig \
    -v
```

In the resulting `output` folder, along with all intermediate results (ROOT files from the GATE simulations) will be two files `eloss_to_wepl_fit_deg3.txt` and `tof_to_wepl_fit_deg3.txt` that contain the coefficients of the polynomials. These files can then be passed to `pctpairprotons` with the corresponding `--fit-kind` parameters.

Here is an example of a `pctpairprotons` invocations with energy-loss fit:
```bash
pctpairprotons \
    -i PhaseSpaceIn_0.root \
    -j PhaseSpaceOut_0.root \
    -o pairs.mhd \
    --plane-in -110 \
    --plane-out 110 \
    --fit output/eloss_to_wepl_fit_deg3.txt \
    --fit-kind energy
```

## Double 1D calibration curve

PCT provides the application `pctdoublelut` that generates two 1D calibration curves from TOF to WEPL. The resulting polynomial fits can then be used when pairing protons with `pctpairprotons` to convert the TOF directly to WEPL. This conversion is not available for energy loss measurements.

Below is an example of how to run `pctdoublelut` to generate calibration curves for a 200 MeV beam:
```bash
pctdoublelut \
    -o output \
    -e 200 \
    -v
```

In the resulting `output` folder, along with all intermediate results (ROOT files from the GATE simulations) will be files named `tof_coeffs_X.txt` and `vel_coeffs_X.txt` that contain the coefficients of the polynomials, where `X` is the polynomial degree. These files can then be passed to `pctpairprotons` with the appropriate parameters.

Here is an example of a `pctpairprotons` invocations that uses the double LUT fits with polynomials of degree 9:
```bash
pctpairprotons \
    -i PhaseSpaceIn_0.root \
    -j PhaseSpaceOut_0.root \
    -o pairs.mhd \
    --plane-in -110 \
    --plane-out 110 \
    --lut-tof output/tof_coeffs_9.txt \
    --lut-vel output/vel_coeffs_9.txt \
    --quadric 1 0 1 0 0 0 0 0 0 -10000 \
    --angle 0
```

Additional notes:
- This method requires the knowledge of the (potentially approximate) support of the object, which is passed to the script as a quadric shape. For additional details about how to define quadric shapes, you can refer to [the documentation of `rtkQuadricShape`](https://www.openrtk.org/Doxygen/classrtk_1_1QuadricShape.html).
- This method can only process list-mode files that contain data for a single acquisition angle, which is specified using the `--angle` parameter. It cannot process data that contain more than one angle.
