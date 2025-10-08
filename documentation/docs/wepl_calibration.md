# WEPL calibration

PCT provides the application `pctweplfit` that generates calibration curves from energy loss and TOF to WEPL. The resulting polynomial fit can then be used when pairing protons with `pctpairprotons` to convert the energy loss or the TOF directly to WEPL. Please note that `pctweplfit` requires GATE, which can be installed following the instructions in the [installation guide](installation.md).

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
    --fit output/tof_to_wepl_fit_deg3.json \
    --fit-kind tof
```

The resulting pairs are directly associated to the corresponding WEPL following the convention explained in the [PCT data format](pct_format.md), that is that $e_\text{in}=0$ and $e_\text{out}=\text{WEPL}$.
