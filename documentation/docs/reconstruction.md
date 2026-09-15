# Reconstruction

PCT supports either distance-driven FDK reconstruction [\[Simon Rit et al, Med. Phys. 2013\]](https://doi.org/10.1118/1.4789589), or gradient descent reconstruction.

## Gradient descent reconstruction

The application that performs gradient descent reconstruction is `pctgradientdescent`. It iteratively solves the reconstruction problem directly from list-mode data. It is much more computationally expensive than distance-driven FDK, and for this reason, runs better on GPU. The GPU part of `pctgradientdescent` was written with the help of [Pytorch](https://pytorch.org/), thus requires the `torch` package to be available in order to run.

A minimal invocation looks like
```bash
pctgradientdescent \
    -p pairs/ \
    -r pairs.*.mhd \
    -g geometry.xml \
    --sp-fit sp_fit.txt
```
where `sp_fit.txt` is a file that maps the proton energy to the corresponding stopping power of water, necessary to compute how much the protons loose energy when traversing the object. This mapping file can be generated using the `pctstoppingpower` application:
```bash
pctstoppingpower -o sp_fit.txt
```

As for all PCT applications, the parameters of `pctgradientdescent` and `pctstoppingpower` can be listed using the `--help` flag:
```bash
pctgradientdescent --help
pctstoppingpower --help
```
