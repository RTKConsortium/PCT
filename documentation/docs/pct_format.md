# PCT list-mode data format

PCT uses its own data format to store list-mode proton CT data as [MetaImage](https://itk.org/Wiki/ITK/MetaIO/Documentation) files (`.mhd` or `.mha` extension). The data corresponds to so-called proton pairs information (position, direction, energy, ...), i.e. the information recorded by the entrance and exit detectors which are paired, e.g.  by a coincidence electronics or by an algorithm.

Some PCT executables produce such data, e.g. `pctpairprotons` from [ROOT](https://root.cern/) data produced by a [GATE](https://github.com/OpenGATE/opengate) simulation, whereas other ones take such data as input, e.g. `pctbinning`. PCT uses an [ITK](https://itk.org/) image internally (`itk::Image`) to read from / write to these [MetaImage](https://itk.org/Wiki/ITK/MetaIO/Documentation) files.

## Data format

A PCT list-mode dataset is a two-dimensional $N\times M$ MetaImage where
* $N$ is the number of proton pairs, i.e. each row corresponds to one proton pair, and
* $M$ is the number of columns in the dataset.

More precisely, the header contains a series of key/value pairs that specify at which column indexes a given information can be found. For instance, if the header specifies that `UpstreamPositionU = 2`, then the third element of each proton pair corresponds to the $u$ coordinate of the proton as detected by the upstream detector.

**Important note:** the data is stored as an ITK image file to facilitate its processing using ITK, but it does **not** represent an image! For this reason, some fields that typically appear in usual header files are meaningless here, such as `ElementSpacing` or `Offset`. Their values are not interpreted in any way.

Below is a complete description of the fields that PCT can handle. For backward-compatibility, some fields have a default value which corresponds to the legacy format (detailed below).

| Header key             | Description                                                  | Default value |
|------------------------|--------------------------------------------------------------|---------------|
| `UpstreamPositionU`    | $u$ coordinate of the proton at the upstream detector        | 0             |
| `UpstreamPositionV`    | $v$ coordinate of the proton at the upstream detector        | 1             |
| `UpstreamPositionW`    | $w$ coordinate of the proton at the upstream detector        | 2             |
| `DownstreamPositionU`  | $u$ coordinate of the proton at the downstream detector      | 3             |
| `DownstreamPositionV`  | $v$ coordinate of the proton at the downstream detector      | 4             |
| `DownstreamPositionW`  | $w$ coordinate of the proton at the downstream detector      | 5             |
| `UpstreamDirectionU`   | Direction along $u$ of the proton at the upstream detector   | 6             |
| `UpstreamDirectionV`   | Direction along $v$ of the proton at the upstream detector   | 7             |
| `UpstreamDirectionW`   | Direction along $w$ of the proton at the upstream detector   | 8             |
| `DownstreamDirectionU` | Direction along $u$ of the proton at the downstream detector | 9             |
| `DownstreamDirectionV` | Direction along $v$ of the proton at the downstream detector | 10            |
| `DownstreamDirectionW` | Direction along $w$ of the proton at the downstream detector | 11            |
| `UpstreamEnergy`       | Kinetic energy of the proton at the upstream detector        | 12            |
| `DownstreamEnergy`     | Kinetic energy of the proton at the downstream detector      | 13            |
| `TrackID`              | `TrackID` of the proton                                      | 14            |
| `WEPL`                 | Water-equivalent path length                                 | ‒             |
| `CreatorProcess`       | Process which created the exit particle                      | ‒             |
| `NuclearProcess`       | Whether the particle encountered a nuclear interaction       | ‒             |
| `Order`                | Number of particle interactions                              | ‒             |
| `TOF`                  | Time-of-flight between the two detectors                     | ‒             |

## Legacy data format

PCT used to use a different data format, which was discontinued due to its rigidity. The new data format should be backward-compatible with files that use the old data format. The description below provides details about this legacy data format.

Proton pairs are stored in 2D images of 3D float vectors, in which the first dimension is a series of 5 or 6 vectors of 3 floats each, one per proton pair, and the second dimension corresponds to the number of proton pairs. Each vector stores the following data:

1. $(u_{\text{in}},v_{\text{in}},w_{\text{in}})$ position of the proton at the entrance detector;
2. $(u_{\text{out}},v_{\text{out}},w_{\text{out}})$ position of the proton at the exit detector;
3. $(\dot u_{\text{in}}, \dot v_{\text{in}}, \dot w_{\text{in}})$ direction of the proton at the entrance detector (as unit vector);
4. $(\dot u_{\text{out}}, \dot v_{\text{out}}, \dot w_{\text{out}})$ direction of the proton at the exit detector (as unit vector);
5. $(e_{\text{in}},e_{\text{out}},t)$ where
    - $e_{\text{in}}$ and $e_{\text{out}}$ are the proton energy at the entrance and exit detectors, respectively.
        - **Important note 1: if $e_{\text{in}}=0$, then $e_{\text{out}}$ will directly be interpreted as the water equivalent path length (WEPL) for the proton pair.**
        - **Important note 2: $e_{\text{in}}$ and $e_{\text{out}}$ can alternatively represent the time of detection of the protons, if `pctpairprotons` was executed using the `--store-time` flag.**
    - $t$ is not used in PCT but can be used to store some useful scalar such as the [GATE](https://github.com/OpenGATE/opengate) `TrackID` in `pctpairprotons`.
6. Optionally, mainly for the work in http://doi.org/10.1088/0031-9155/61/9/3258 $(\text{creatorProcess}, \text{nuclearProcess}, \text{order})$ are used in simulations to indicate
    - the process which created the exit particle,
    - whether the particle did encounter a nuclear interaction,
    - the number of particle interactions.

Note that PCT assumes that the proton beam goes along the $w$ axis in the positive direction ($w_{\text{in}} < w_{\text{out}}$).
