#!/usr/bin/env python
import argparse
from datetime import datetime
import time
import struct
import os
import sys
import random
import math
import itk
from itk import PCT as pct
from tqdm import tqdm
import torch
from numpy import loadtxt, flip

m0 = 938.27208943  # MeV/c2
c = 299.792458  # mm/ns

LMComponentType = itk.ctype("float")
LMPixelType = itk.Vector[LMComponentType, 3]
LMImageType = itk.Image[LMPixelType, 2]

RSPPixelType = itk.ctype("double")
RSPImageType = itk.Image[RSPPixelType, 3]


def build_parser():
    parser = pct.PCTArgumentParser(
        description="Gradient descent reconstruction of PCT data",
    )
    parser.add_argument(
        "--path", "-p", help="Path containing list-mode files", type=str, required=True
    )
    parser.add_argument(
        "--regexp",
        "-r",
        help="Regular expression to select list-mode files in path",
        type=str,
        required=True,
    )
    parser.add_argument("-g", "--geometry", help="Geometry file", required=True)
    parser.add_argument("-o", "--output-dir", help="Output directory")
    parser.add_argument(
        "--sp-fit", help="Fit of stopping power of water", required=True
    )
    parser.add_argument("--rsp-map", help="Initial relative stopping power map")
    parser.add_argument(
        "--initial-energy",
        help="Initial energy of the protons in MeV",
        type=float,
        default=200.0,
    )
    parser.add_argument(
        "--step", help="Step size (in mm along z)", type=float, default=1.0
    )
    parser.add_argument(
        "--polydeg", help="Polynomial degree for MLP estimation", type=int, default=5
    )
    parser.add_argument(
        "--hull",
        "--quadric",
        help="Quadric representing the hull of the object",
        type=float,
        nargs=10,
        required=False,
    )
    parser.add_argument("--optimizer", help="Optimizer", default="Adam")
    parser.add_argument(
        "-n", "--number-of-iterations", help="Number of iterations", type=int, default=5
    )
    parser.add_argument(
        "-l", "--learning-rate", help="Learning rate", type=float, default=0.3
    )
    parser.add_argument(
        "--tv", help="Total variation regularization", type=float, default=0.0
    )
    parser.add_argument(
        "--max-batch-size", help="Batch size", type=int, required=False, default=100_000
    )
    parser.add_argument(
        "-w",
        "--num-workers",
        help="Number of workers to compute MLPs",
        type=int,
        required=False,
        default=0,
    )
    parser.add_argument(
        "-q",
        "--physical-quantity",
        help="Physical quantity to use to do the reconstruction",
        default="time",
        choices=["energy", "time"],
    )
    parser.add_argument(
        "--cpu",
        help="Use CPU even if a GPU is available",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-s",
        "--number-of-subsets",
        help="Number of batches in a subset",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--verbose", "-v", help="Verbose execution", default=False, action="store_true"
    )

    pct3Doutputimage_group = parser.add_argument_group("Output 3D image properties")
    pct3Doutputimage_group.add_argument(
        "--size", help="Size", type=int, nargs="+", default=[110, 3, 110]
    )
    pct3Doutputimage_group.add_argument(
        "--spacing", help="Spacing", type=float, nargs="+", default=[2.0, 10.0, 2.0]
    )
    return parser


def pv(verbosity, msg):
    if verbosity:
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")


class ProtonDataset(torch.utils.data.Dataset):

    def __init__(self, inputs, angles, polydeg, step, hull, pq):

        self.inputs = inputs
        self.angles = angles
        self.polydeg = polydeg
        self.step = step
        self.hull = hull
        self.pq = pq

        idx_bins = [0]
        raw_files = []
        for i in inputs:
            # Read image size without loading it
            reader = itk.MetaImageIO.New()
            reader.SetFileName(i)
            reader.ReadImageInformation()
            number_of_elements = reader.GetDimensions(1)
            idx_bins.append(idx_bins[-1] + number_of_elements)
            raw_files.append(ProtonDataset.get_mhd_value(i, "ElementDataFile"))

        self.idx_bins = torch.tensor(idx_bins)
        self.raw_files = raw_files

        self.element_size = 5 * 3  # bytes

        mhd_to_bytes = {
            "MET_CHAR": 1,
            "MET_UCHAR": 1,
            "MET_SHORT": 2,
            "MET_USHORT": 2,
            "MET_INT": 4,
            "MET_UINT": 4,
            "MET_LONG_LONG": 8,
            "MET_ULONG_LONG": 8,
            "MET_FLOAT": 4,
            "MET_DOUBLE": 8,
        }
        element_type = ProtonDataset.get_mhd_value(
            inputs[0], "ElementType"
        )  # assuming they are all the same across files
        self.pixel_size = mhd_to_bytes[element_type]

    @staticmethod
    def get_mhd_value(mhd_filename, requested_key):
        with open(mhd_filename) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                key, value = line.split("=", 1)
                if key.strip().lower() == requested_key.lower():
                    return value.strip()
        return None

    def get_idx_file(self, idx):
        return (torch.bucketize(idx, self.idx_bins, right=True) - 1).item()

    def get_element(self, idx):
        idx_file = self.get_idx_file(idx)
        idx_element = idx - self.idx_bins[idx_file]
        raw_file = os.path.join(
            os.path.dirname(self.inputs[idx_file]), self.raw_files[idx_file]
        )
        with open(raw_file, "rb") as f:
            f.seek(idx_element * self.element_size * self.pixel_size)
            r = f.read(self.element_size * self.pixel_size)
            return struct.unpack("f" * self.element_size, r)

    def __len__(self):
        return self.idx_bins[-1]

    def __getitem__(self, idx):

        element = self.get_element(idx)

        pos_in = [element[0], element[1], element[2]]
        pos_out = [element[3], element[4], element[5]]
        dir_in = [element[6], element[7], element[8]]
        dir_out = [element[9], element[10], element[11]]

        if self.pq == "time":
            pq = element[13] - element[12]
        elif self.pq == "energy":
            pq = element[13]
        else:
            raise ValueError

        ws = torch.arange(pos_in[2], pos_out[2] + self.step, self.step)

        if self.hull is None:

            mlp = pct.PolynomialMLPFunction.New()
            mlp.SetPolynomialDegree(self.polydeg)
            mlp.Init(pos_in, pos_out, dir_in, dir_out)
            us, vs = mlp.Evaluate(ws.tolist())
            us, vs = torch.tensor(us), torch.tensor(vs)

        else:

            hull = rtk.QuadricShape.New()
            hull.SetA(self.hull[0])
            hull.SetB(self.hull[1])
            hull.SetC(self.hull[2])
            hull.SetD(self.hull[3])
            hull.SetE(self.hull[4])
            hull.SetF(self.hull[5])
            hull.SetG(self.hull[6])
            hull.SetH(self.hull[7])
            hull.SetI(self.hull[8])
            hull.SetJ(self.hull[9])

            dir_out_r = [-dir_out[0], -dir_out[1], -dir_out[2]]
            intersect_u, d_in, _ = hull.IsIntersectedByRay(pos_in, dir_in)
            intersect_d, d_out, _ = hull.IsIntersectedByRay(pos_out, dir_out_r)
            if not (intersect_u and intersect_d):
                us = torch.linspace(pos_in[0], pos_out[0], ws.shape[0])
                vs = torch.linspace(pos_in[1], pos_out[1], ws.shape[0])
            else:

                int_in = [pos_in[i] + dir_in[i] * d_in for i in range(3)]
                int_out = [pos_out[i] - dir_out[i] * d_out for i in range(3)]

                us = torch.empty_like(ws)
                vs = torch.empty_like(ws)

                # Before object
                mask_before = ws < int_in[2]
                d_before = (ws[mask_before] - pos_in[2]) / (int_in[2] - pos_in[2])
                us[mask_before] = pos_in[0] + dir_in[0] * d_in * d_before
                vs[mask_before] = pos_in[1] + dir_in[1] * d_in * d_before

                # After object
                mask_after = ws > int_out[2]
                d_after = (ws[mask_after] - pos_out[2]) / (int_out[2] - pos_out[2])
                us[mask_after] = pos_out[0] - dir_out[0] * d_out * d_after
                vs[mask_after] = pos_out[1] - dir_out[1] * d_out * d_after

                # In object
                mask_in = torch.logical_not(torch.logical_or(mask_before, mask_after))
                mlp = pct.PolynomialMLPFunction.New()
                mlp.SetPolynomialDegree(self.polydeg)
                mlp.Init(int_in, int_out, dir_in, dir_out)
                us_in, vs_in = mlp.Evaluate(ws[mask_in].tolist())
                us[mask_in], vs[mask_in] = torch.tensor(
                    us_in, device=us.device
                ), torch.tensor(vs_in, device=vs.device)

        angle = self.angles[self.get_idx_file(idx)]

        return pq, angle, (us, vs, ws)


class RandomSampler:

    def __init__(self, length):
        self.length = length

    def __len__(self):
        return self.length

    def __iter__(self):
        yield from self.random_range(0, self.length, 1)

    @staticmethod
    def random_range(start, stop=None, step=None):
        # https://stackoverflow.com/a/53551417
        if stop == None:
            start, stop = 0, start
        if step == None:
            step = 1
        mapping = lambda i: (i * step) + start
        maximum = (stop - start) // step
        value = random.randint(0, maximum)
        offset = random.randint(0, maximum) * 2 + 1
        multiplier = 4 * (maximum // 4) + 1
        modulus = int(2 ** math.ceil(math.log2(maximum)))
        found = 0
        while found < maximum:
            if value < maximum:
                found += 1
                yield mapping(value)
            value = (value * multiplier + offset) % modulus


class Velocity(torch.autograd.Function):

    @staticmethod
    def forward(ctx, e):
        ctx.save_for_backward(e)
        return c * torch.sqrt(1 - (m0**2) / ((m0 + e) ** 2))

    @staticmethod
    def backward(ctx, grad_output):
        (e,) = ctx.saved_tensors
        dv_de = c * m0**2 / ((e + m0) ** 3 * torch.sqrt(-(m0**2) / (e + m0) ** 2 + 1))
        return grad_output * dv_de


class EnergyLoss(torch.autograd.Function):

    @staticmethod
    def forward(ctx, es, rsps, sp_h2o, d):
        ctx.save_for_backward(sp_h2o, d)
        sp = rsps * sp_h2o
        eloss = sp * d
        es2 = es - eloss
        return es2

    @staticmethod
    def backward(ctx, grad_output):
        sp_h2o, d = ctx.saved_tensors
        grad_rsps = -sp_h2o * d * grad_output
        return None, grad_rsps, None, None


class RSPInterp(torch.autograd.Function):

    @staticmethod
    def forward(ctx, rsp, rsp_spacing, x, y, z):
        upper_bounds = [
            (rsp.shape[i] * rsp_spacing[i] / 2) - rsp_spacing[i] / 2 for i in range(3)
        ]
        lower_bounds = [-upper_bounds[i] for i in range(3)]

        xs = torch.linspace(
            lower_bounds[0], upper_bounds[0], rsp.shape[0], device=x.device
        )
        ys = torch.linspace(
            lower_bounds[1], upper_bounds[1], rsp.shape[1], device=y.device
        )
        zs = torch.linspace(
            lower_bounds[2], upper_bounds[2], rsp.shape[2], device=z.device
        )

        x_idx = torch.searchsorted(xs, x)
        y_idx = torch.searchsorted(ys, y)
        z_idx = torch.searchsorted(zs, z)

        in_mask = (
            (0 < x_idx)
            & (x_idx < xs.shape[0])
            & (0 < y_idx)
            & (y_idx < ys.shape[0])
            & (0 < z_idx)
            & (z_idx < zs.shape[0])
        )

        oob_value = 0.0
        cs = torch.full((len(x),), oob_value, device=x.device)

        x_in = x[in_mask]
        y_in = y[in_mask]
        z_in = z[in_mask]

        x_idx_in = x_idx[in_mask]
        y_idx_in = y_idx[in_mask]
        z_idx_in = z_idx[in_mask]

        x0, x1 = xs[x_idx_in - 1], xs[x_idx_in]
        y0, y1 = ys[y_idx_in - 1], ys[y_idx_in]
        z0, z1 = zs[z_idx_in - 1], zs[z_idx_in]
        xd = (x_in - x0) / (x1 - x0)
        yd = (y_in - y0) / (y1 - y0)
        zd = (z_in - z0) / (z1 - z0)
        c000 = rsp[x_idx_in - 1, y_idx_in - 1, z_idx_in - 1]
        c001 = rsp[x_idx_in - 1, y_idx_in - 1, z_idx_in]
        c010 = rsp[x_idx_in - 1, y_idx_in, z_idx_in - 1]
        c011 = rsp[x_idx_in - 1, y_idx_in, z_idx_in]
        c100 = rsp[x_idx_in, y_idx_in - 1, z_idx_in - 1]
        c101 = rsp[x_idx_in, y_idx_in - 1, z_idx_in]
        c110 = rsp[x_idx_in, y_idx_in, z_idx_in - 1]
        c111 = rsp[x_idx_in, y_idx_in, z_idx_in]
        c00 = c000 * (1 - xd) + c100 * xd
        c01 = c001 * (1 - xd) + c101 * xd
        c10 = c010 * (1 - xd) + c110 * xd
        c11 = c011 * (1 - xd) + c111 * xd
        c0 = c00 * (1 - yd) + c10 * yd
        c1 = c01 * (1 - yd) + c11 * yd

        cs[in_mask] = c0 * (1 - zd) + c1 * zd

        ctx.rsp_shape = rsp.shape
        ctx.save_for_backward(in_mask, xd, yd, zd, x_idx_in, y_idx_in, z_idx_in)

        return cs

    @staticmethod
    def backward(ctx, grad_output):
        in_mask, xd, yd, zd, x_idx_in, y_idx_in, z_idx_in = ctx.saved_tensors
        rsp_shape = ctx.rsp_shape

        grad_output_in = grad_output[in_mask]
        drsp = torch.zeros(rsp_shape, device=grad_output.device)

        for xx in [False, True]:
            for yy in [False, True]:
                for zz in [False, True]:

                    x = x_idx_in if xx else x_idx_in - 1
                    y = y_idx_in if yy else y_idx_in - 1
                    z = z_idx_in if zz else z_idx_in - 1
                    idx_lin = z + y * rsp_shape[2] + x * rsp_shape[2] * rsp_shape[1]

                    xd_ = xd if xx else 1 - xd
                    yd_ = yd if yy else 1 - yd
                    zd_ = zd if zz else 1 - zd

                    drsp.put_(
                        idx_lin, grad_output_in * xd_ * yd_ * zd_, accumulate=True
                    )

        return drsp, None, None, None, None


class SPWater(torch.autograd.Function):

    @staticmethod
    def forward(ctx, energies, sps, es):

        es_idx = torch.searchsorted(energies, es)

        # Deal with incorrect values due to wrong energies
        # Just to avoid errors, as values are discarded later on anyway
        es_idx[es_idx <= 0] = 0
        es_idx[es_idx >= energies.shape[0]] = energies.shape[0] - 1

        sp0, sp1 = sps[es_idx - 1], sps[es_idx]
        e0, e1 = energies[es_idx - 1], energies[es_idx]
        ed = (e1 - es) / (e1 - e0)
        cs = ed * sp0 + (1 - ed) * sp1

        ctx.save_for_backward(e0, e1, sp0, sp1)

        return cs

    @staticmethod
    def backward(ctx, grad_output):
        e0, e1, sp0, sp1 = ctx.saved_tensors
        des = (sp1 / (e1 - e0)) - (sp0 / (e1 - e0))
        return None, None, des * grad_output


def forward_iter(w, xs, ys, zs, es, rsp, rsp_spacing, spfit_data):

    x, y, z = xs[:, w].contiguous(), ys[:, w].contiguous(), zs[:, w].contiguous()
    x2, y2, z2 = xs[:, w + 1], ys[:, w + 1], zs[:, w + 1]
    d = torch.sqrt((x - x2) ** 2 + (y - y2) ** 2 + (z - z2) ** 2)

    rsps = RSPInterp.apply(rsp, rsp_spacing, x, y, z)
    sp_h2o = SPWater.apply(
        spfit_data[:, 0].contiguous(), spfit_data[:, 1].contiguous() / 10.0, es
    )
    sp = rsps * sp_h2o
    eloss = sp * d
    es2 = es - eloss

    # Avoid NaNs during execution, but protons need to be filtered out afterwards
    es2 = torch.where(es2 < 1.0, 1.0, es2)

    vel = Velocity.apply(es)
    vel2 = Velocity.apply(es2)

    # Trapezoidal rule
    tof = (2.0 / (vel + vel2)) * d

    return tof, es2


def forward_model(rsp, rsp_spacing, mlps, angles, spfit_data, e0):

    us, vs, ws = mlps
    assert us.shape == vs.shape

    # Apply rotation
    thetas = angles.unsqueeze(1)
    xs = us * torch.cos(thetas) + ws * torch.sin(thetas)
    ys = vs
    zs = -us * torch.sin(thetas) + ws * torch.cos(thetas)

    number_of_protons = xs.shape[0]
    es = torch.full((number_of_protons,), e0, device=us.device)
    tofs = torch.full((number_of_protons,), 0.0, device=us.device)

    for w_idx, _ in enumerate(ws[0, :-1]):
        tof, es2 = forward_iter(w_idx, xs, ys, zs, es, rsp, rsp_spacing, spfit_data)

        tofs += tof
        es = es2

    return tofs, es


def total_variation(img):
    dx = torch.sum(torch.abs(img[1:, :, :] - img[:-1, :, :]))
    dy = torch.sum(torch.abs(img[:, 1:, :] - img[:, :-1, :]))
    dz = torch.sum(torch.abs(img[:, :, 1:] - img[:, :, :-1]))
    return dx + dy + dz


def save_img(img, spacing, filename):
    origin = tuple(-(img.shape[i] * spacing[i] / 2) + spacing[i] / 2 for i in range(3))

    img_itk = itk.GetImageFromArray(img)
    img_itk.SetSpacing(spacing)
    img_itk.SetOrigin(origin)
    itk.imwrite(img_itk, filename)


def write_report(
    output_dir,
    iteration="iteration",
    subset="subset",
    loss="loss",
    subset_time="subset_time",
    backward_time="backward_time",
    iteration_elapsed="iteration_elapsed",
    total_elapsed="total_elapsed",
):
    with open(os.path.join(output_dir, "report.txt"), "a", encoding="utf-8") as f:
        f.write(
            "\t".join(
                map(
                    str,
                    [
                        iteration,
                        subset,
                        loss,
                        subset_time,
                        backward_time,
                        iteration_elapsed,
                        total_elapsed,
                    ],
                )
            )
            + os.linesep
        )


def get_optimizer(name, learning_rate, rsp):
    if name == "SGD":
        optimizer = torch.optim.SGD([rsp], lr=learning_rate)
    elif name == "Adagrad":
        optimizer = torch.optim.Adagrad([rsp], lr=learning_rate)
    elif name == "Adam":
        optimizer = torch.optim.Adam([rsp], lr=learning_rate)
    else:
        raise NotImplementedError
    optimizer.zero_grad()
    return optimizer


def process(args_info: argparse.Namespace):

    from itk import RTK as rtk

    torch.multiprocessing.set_start_method("spawn", force=True)
    torch.multiprocessing.set_sharing_strategy(
        "file_system"
    )  # https://stackoverflow.com/a/73289157
    torch.autograd.set_detect_anomaly(True, check_nan=True)

    device = torch.device(
        "cuda:0" if torch.cuda.is_available() and not args_info.cpu else "cpu"
    )
    pv(args_info.verbose, f"Will use device {device} for computations")

    spfit_data = torch.tensor(flip(loadtxt(args_info.sp_fit, skiprows=4), 0).copy()).to(
        device
    )
    pv(args_info.verbose, "Read SP fit data.")

    pv(args_info.verbose, "Attempt at reading checkpoint…")
    os.makedirs(args_info.output_dir, exist_ok=True)
    checkpoint_path = os.path.join(args_info.output_dir, "checkpoint.pth")
    try:
        checkpoint = torch.load(checkpoint_path, weights_only=False)
        iteration = checkpoint["iteration"]
        rsp = checkpoint["rsp"].to(device)
        rsp_spacing = checkpoint["rsp_spacing"]
        optimizer = checkpoint["optimizer"]
        loss = checkpoint["loss"]
        rsp_shape = rsp.shape
        pv(args_info.verbose, "Checkpoint successfully read!")
    except FileNotFoundError:
        pv(args_info.verbose, "No checkpoint found, starting from scratch.")

        if args_info.rsp_map is None:
            rsp_shape = args_info.size
            rsp_spacing = args_info.spacing
            rsp = torch.full(rsp_shape, 0.0, device=device)
        else:
            rsp_itk = itk.imread(rsp_map)
            rsp = torch.from_numpy(itk.GetArrayFromImage(rsp_itk)).to(device)
            rsp_shape = rsp.shape
            rsp_spacing = tuple(rsp_itk.GetSpacing())
        pv(
            args_info.verbose,
            f"Initialized RSP map with shape {rsp_shape} and spacing {rsp_spacing}.",
        )

        optimizer = get_optimizer(args_info.optimizer, args_info.learning_rate, rsp)
        iteration = 1
        loss = torch.nn.MSELoss(reduction="mean")
        write_report(args_info.output_dir)
    loss_tot = []

    rsp.requires_grad = True

    pv(args_info.verbose, "Reading inputs and geometry…")
    names = itk.RegularExpressionSeriesFileNames.New()
    names.SetDirectory(args_info.path)
    names.SetNumericSort(False)
    names.SetRegularExpression(args_info.regexp)
    names.SetSubMatch(0)
    inputs = names.GetFileNames()
    pv(args_info.verbose, f"Regular expression matches {len(inputs)} file(s)…")

    pv(args_info.verbose, "Reading geometry…")
    geometryReader = rtk.ThreeDCircularProjectionGeometryXMLFileReader.New()
    geometryReader.SetFilename(args_info.geometry)
    geometryReader.GenerateOutputInformation()
    angles = torch.tensor(geometryReader.GetGeometry().GetSourceAngles())
    angles -= torch.pi / 2.0  # realign with rest of PCT
    nangles = len(angles)
    ninputs = len(inputs)
    assert nangles == ninputs, f"{nangles} angles, but {ninputs} inputs!"

    pv(args_info.verbose, "Creating dataset…")
    dataset = ProtonDataset(
        inputs,
        angles,
        args_info.polydeg,
        args_info.step,
        args_info.hull,
        args_info.physical_quantity,
    )
    number_of_protons = len(dataset)
    protons_per_subset = number_of_protons // args_info.number_of_subsets
    # Calculate "ideal" batch size: largest divisor below some threshold
    batch_size = protons_per_subset
    batch_per_subset = 1
    while batch_size > args_info.max_batch_size:
        batch_per_subset += 1
        batch_size = protons_per_subset // batch_per_subset
    dataloader = torch.utils.data.DataLoader(
        dataset,
        sampler=RandomSampler(len(dataset)),
        batch_size=batch_size,
        num_workers=args_info.num_workers,
        in_order=False,
        persistent_workers=args_info.num_workers > 0,
        pin_memory=True,
        drop_last=False,
    )
    pv(
        args_info.verbose,
        f"{number_of_protons} protons split into {args_info.number_of_subsets} subsets of {protons_per_subset} protons with {batch_per_subset} batches of {batch_size} protons per subset.",
    )

    start_time = time.time()

    while iteration <= args_info.number_of_iterations:
        pv(args_info.verbose, f"Starting iteration {iteration}…")

        dataloader_it = iter(dataloader)

        iteration_start_time = time.time()

        for subset in tqdm(range(1, args_info.number_of_subsets + 1)):

            subset_start_time = time.time()
            backward_time = 0

            optimizer.zero_grad()
            loss_tot = []

            pq_cpu, angle_cpu, mlp_cpu = next(dataloader_it)

            for _ in range(batch_per_subset):

                angle = angle_cpu.to(device)
                mlp = (
                    mlp_cpu[0].to(device),
                    mlp_cpu[1].to(device),
                    mlp_cpu[2].to(device),
                )

                tofs_fm, es_fm = forward_model(
                    rsp, rsp_spacing, mlp, angle, spfit_data, args_info.initial_energy
                )

                if args_info.physical_quantity == "time":
                    pq_fm = tofs_fm
                elif args_info.physical_quantity == "energy":
                    pq_fm = es_fm
                else:
                    raise ValueError

                pq_m = pq_cpu.to(device, dtype=pq_fm.dtype)

                proton_mask = es_fm > 1.0
                loss_n = loss(
                    pq_fm[proton_mask], pq_m[proton_mask]
                ) + args_info.tv * total_variation(rsp)
                loss_tot.append(loss_n.item())

                backward_start_time = time.time()
                loss_n.backward()
                backward_time += time.time() - backward_start_time

            optimizer.step()

            # Enforce positivity in RSP
            with torch.no_grad():
                rsp[rsp < 0.0] = 0.0

            loss_mean = torch.tensor(loss_tot).mean().item()
            subset_time = time.time() - subset_start_time
            iteration_elapsed = time.time() - iteration_start_time
            total_elapsed = time.time() - start_time
            write_report(
                args_info.output_dir,
                iteration=iteration,
                subset=subset,
                loss=loss_mean,
                subset_time=subset_time,
                backward_time=backward_time,
                iteration_elapsed=iteration_elapsed,
                total_elapsed=total_elapsed,
            )
            save_img(
                rsp.cpu().detach().numpy(),
                rsp_spacing,
                os.path.join(
                    args_info.output_dir, f"iteration_{iteration}_subset_{subset}.mhd"
                ),
            )

        iteration += 1

        checkpoint = {
            "iteration": iteration,
            "rsp": rsp,
            "rsp_spacing": rsp_spacing,
            "optimizer": optimizer,
            "loss": loss,
        }
        torch.save(checkpoint, checkpoint_path)

    end_time = time.time()

    pv(
        args_info.verbose,
        f"{args_info.number_of_iterations} iterations completed in {end_time - start_time} s.",
    )


def main(argv=None):
    parser = build_parser()
    args_info = parser.parse_args(argv)
    process(args_info)


if __name__ == "__main__":
    main()
