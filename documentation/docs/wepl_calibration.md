# WEPL calibration

PCT provides two applications that generate calibration curves to convert energy-loss or TOF to the corresponding WEPL. These calibration curves are based on the works of [Ulrich-Pur et al. (2023)](https://doi.org/10.1088/1748-0221/18/02/C02062) (`pctsdpweplfit`) and [Coussat et al. (2025)](https://hal.science/hal-05375602) (`pctweplfit`), respectively.

## `pctsdpweplfit`

Based on the method described by [Ulrich-Pur et al. (2023)](doi.org/10.1088/1748-0221/18/02/C02062), this application creates a calibration curve for slowing-down power to WEPL conversion. The resulting polynomial fit can then be used when pairing protons with `pctpairprotons` to convert the ΔTOF to WEPL. This method is not compatible with energy-loss measurements.

Below is an example of how to run `pctsdpweplfit`:
```bash
pctsdpweplfit \
    -o output \
    --wepl-samples 20 \
    --max-wepl 200 \
    --number-of-particles 100 \
    --detector-distance 350 \
    --initial-energy 200 \
    -v
```

In the resulting `output` folder, along with all intermediate results (ROOT files from the GATE simulations) will be the file `fit.json` that contains the coefficients of the polynomial.

## `pctweplfit`

Based on the method described by [Coussat et al. (2025)](https://hal.science/hal-05375602), this application generates calibration curves from energy loss and TOF to WEPL. The resulting polynomial fit can then be used when pairing protons with `pctpairprotons` to convert the energy loss or the TOF directly to WEPL.

Below is an example of how to run `pctweplfit`:
```bash
pctweplfit \
    -o output \
    --savefig \
    -v
```

In the resuling `output` folder, along with all intermediate results (ROOT files from the GATE simulations) will be two files `eloss_to_wepl_fit_deg3.json` and `tof_to_wepl_fit_deg3.json` that contain the coefficients for the polynomials. These files can then be passed to `pctpairprotons` with the corresponding `--fit-kind` parameters.

Here is an example of a `pctpairprotons` invocations with energy-loss fit:
```bash
pctpairprotons \
    -i PhaseSpaceIn_0.root \
    -j PhaseSpaceOut_0.root \
    -o pairs.mhd \
    --plane-in -110 \
    --plane-out 110 \
    --fit output/eloss_to_wepl_fit_deg3.json \
    --fit-kind energy
```

## Additional notes

- The resulting pairs are directly associated to the corresponding WEPL following the convention explained in the [PCT data format](pct_format.md), that is that $e_\text{in}=0$ and $e_\text{out}=\text{WEPL}$.
- These applications require GATE, which can be installed following the instructions in the [installation guide](installation.md).
