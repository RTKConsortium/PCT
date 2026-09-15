#!/usr/bin/env python
import argparse
import json
import itk
from itk import PCT as pct
import numpy as np
import numpy.lib.recfunctions as rfn


def build_parser():
    parser = pct.PCTArgumentParser(
        description="Pair corresponding protons from GATE ROOT files"
    )
    parser.add_argument(
        "-i",
        "--input-in",
        help="Root phase space file of particles before object",
        required=True,
    )
    parser.add_argument(
        "-j",
        "--input-out",
        help="Root phase space file of particles after object",
        required=True,
    )
    parser.add_argument("-o", "--output", help="Output file name", required=True)
    parser.add_argument(
        "--plane-in",
        help="Plane position of incoming protons",
        required=True,
        type=float,
    )
    parser.add_argument(
        "--plane-out",
        help="Plane position of outgoing protons",
        required=True,
        type=float,
    )
    parser.add_argument(
        "--min-run", help="Minimum run (inclusive)", default=0, type=int
    )
    parser.add_argument(
        "--max-run", help="Maximum run (exclusive)", default=1e6, type=int
    )
    parser.add_argument(
        "--no-nuclear",
        help="Remove inelastic nuclear collisions",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--store-time",
        help="Store time instead of energy in the output list-mode",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--verbose", "-v", help="Verbose execution", default=False, action="store_true"
    )
    parser.add_argument(
        "--psin", help="Name of tree in input phase space", default="PhaseSpace"
    )
    parser.add_argument(
        "--psout", help="Name of tree in output phase space", default="PhaseSpace"
    )

    single_lut = parser.add_argument_group("Single LUT for conversion to WEPL")
    single_lut.add_argument(
        "--fit", help="Fit file used to convert from energy loss or TOF to WEPL"
    )
    single_lut.add_argument(
        "--fit-kind",
        help="Whether to convert to WEPL using energy loss or TOF",
        choices=["tof", "energy"],
    )

    double_lut = parser.add_argument_group("Double LUT for conversion to WEPL")
    double_lut.add_argument(
        "--lut-tof", help="LUT corresponding to TOF in double LUT method"
    )
    double_lut.add_argument(
        "--lut-vel", help="LUT corresponding to velocity in double LUT method"
    )
    double_lut.add_argument(
        "--quadric",
        help="Quadric representing the hull of the object",
        type=float,
        nargs=10,
        required=False,
    )
    double_lut.add_argument(
        "--angle", help="Angle of the current projection", type=float, required=False
    )

    return parser


m0 = 938.27208943  # MeV/c2
c = 299.792458  # mm/ns


def make_hull(quadric, angle):
    from itk import RTK as rtk

    hull = rtk.QuadricShape.New()
    hull.SetA(quadric[0])
    hull.SetB(quadric[1])
    hull.SetC(quadric[2])
    hull.SetD(quadric[3])
    hull.SetE(quadric[4])
    hull.SetF(quadric[5])
    hull.SetG(quadric[6])
    hull.SetH(quadric[7])
    hull.SetI(quadric[8])
    hull.SetJ(quadric[9])

    theta = np.deg2rad(angle)
    hull.Rotate(
        itk.matrix_from_array(
            [
                [np.cos(theta), 0.0, np.sin(theta)],
                [0.0, 1.0, 0.0],
                [-np.sin(theta), 0.0, np.cos(theta)],
            ]
        )
    )

    return hull


def tof_in_vacuum(d, e):
    v = c * np.sqrt(1 - (m0**2) / ((m0 + e) ** 2))
    return d / v


def get_tof_and_distance(hull, p_u, d_u, p_d, d_d, energy, tof):

    d_d_r = [-d for d in d_d]

    intersect_u, d_uo, _ = hull.IsIntersectedByRay(p_u, d_u)  # assuming d_u is unitary
    intersect_d, d2, _ = hull.IsIntersectedByRay(p_d, d_d_r)
    if not (intersect_u and intersect_d):
        raise ValueError

    int1 = [p_u[i] + d_u[i] * d_uo for i in range(3)]
    int2 = [p_d[i] + d_d_r[i] * d2 for i in range(3)]
    d1 = np.sqrt(np.sum([(int1[i] - int2[i]) ** 2 for i in range(3)]))

    tof_od = tof - tof_in_vacuum(d_uo, energy)

    return tof_od, d1, d2


def convert_tof_to_wepl(fit_tof, fit_vel, pairs, quadric, angle):
    tof_coeffs = np.loadtxt(fit_tof)
    vel_coeffs = np.loadtxt(fit_vel)

    def tof_to_wepl(tof, d1, d2):
        t1 = d1 * np.polymul(tof_coeffs, vel_coeffs)
        t2 = [d2, 0.0]
        t3 = tof * np.polymul(vel_coeffs, [1.0, 0.0])
        p = np.polysub(np.polyadd(t1, t2), t3)
        roots = [
            np.real(r)
            for r in np.roots(p)
            if np.imag(r) == 0.0
            and np.real(r) >= 0.0
            and (d1 / r) * np.polyval(tof_coeffs, r) < tof
        ]
        return roots[-1]

    names = pairs.dtype.names
    hull = make_hull(quadric, angle)

    for pi, pair in enumerate(pairs):
        p = dict(zip(names, pair.tolist()))
        p_u = [p["u_in"], p["v_in"], p["w_in"]]
        p_d = [p["u_out"], p["v_out"], p["w_out"]]
        d_u = [p["du_in"], p["dv_in"], p["dw_in"]]
        d_d = [p["du_out"], p["dv_out"], p["dw_out"]]
        energy = p["KineticEnergy_in"]
        tof = p["LocalTime_out"] - p["LocalTime_in"]
        try:
            tof, d1, d2 = get_tof_and_distance(hull, p_u, d_u, p_d, d_d, energy, tof)
            wepl = tof_to_wepl(tof, d1, d2)
        except ValueError:
            wepl = 0.0

        pairs[pi]["KineticEnergy_out"] = wepl
        pairs[pi]["KineticEnergy_in"] = 0.0


def process(args_info: argparse.Namespace):
    import uproot

    if args_info.verbose:

        def verbose(message):
            print(message)

    else:

        def verbose(message):
            pass

    measurement_column = "GlobalTime" if args_info.store_time else "KineticEnergy"

    def load_tree_as_df(root_file, tree_name):

        tree = uproot.open(root_file)[tree_name]
        branches = tree.arrays(library="np")

        # Some versions of uproot return a dictionnary, and some others a numpy.ndarray
        # We handle both cases here to generate the dtype
        if isinstance(branches, dict):
            dtype = [(name, branch.dtype) for name, branch in branches.items()]
        elif isinstance(branches, np.ndarray):
            dtype = branches.dtype.descr
        else:
            raise NotImplementedError

        ps = np.rec.recarray((len(branches["RunID"]),), dtype=dtype)
        for branch_name, _ in dtype:
            ps[branch_name] = branches[branch_name]

        ps = rfn.rename_fields(
            ps,
            {
                "Position_X": "u",
                "Position_Y": "v",
                "Position_Z": "w",
            },
        )
        ps = rfn.rename_fields(
            ps,
            {
                "Direction_X": "du",
                "Direction_Y": "dv",
                "Direction_Z": "dw",
            },
        )

        ps = ps[(ps["RunID"] >= args_info.min_run) & (ps["RunID"] < args_info.max_run)]

        return ps

    ps_in = load_tree_as_df(args_info.input_in, args_info.psin)
    ps_in["w"] = args_info.plane_in

    verbose("Read input phase space:\n" + str(ps_in))
    ps_out = load_tree_as_df(args_info.input_out, args_info.psout)
    ps_out["w"] = args_info.plane_out
    verbose("Read output phase space:\n" + str(ps_out))

    merge_columns = ["RunID", "EventID"]
    if args_info.no_nuclear:
        merge_columns.append("TrackID")

    # Remove duplicates
    _, unique_index = np.unique(ps_in[merge_columns], return_index=True)
    ps_in = ps_in[unique_index]

    ps_in.dtype.names = [
        n if n in merge_columns else n + "_in" for n in ps_in.dtype.names
    ]
    ps_out.dtype.names = [
        n if n in merge_columns else n + "_out" for n in ps_out.dtype.names
    ]
    ps_in_uniques = [n for n in ps_in.dtype.names if n not in merge_columns]
    ps_out_uniques = [n for n in ps_out.dtype.names if n not in merge_columns]

    if args_info.no_nuclear:
        # Easy case, there should be at most one row in ps_in and ps_out for keys ['RunID', 'EventID', 'TrackID']
        intersect, intersect_in, intersect_out = np.intersect1d(
            ps_in[merge_columns], ps_out[merge_columns], return_indices=True
        )
        pairs = rfn.merge_arrays(
            (
                intersect,
                ps_in[ps_in_uniques][intersect_in],
                ps_out[ps_out_uniques][intersect_out],
            ),
            asrecarray=True,
            flatten=True,
        )
    else:
        # More complicated, there can be more than one row in ps_out for keys ['RunID', 'EventID'] (because of 'TrackID')
        # The solution is to repeat the computation for many TrackIDs, then merge it all together
        track_max = ps_out["TrackID_out"].max()
        verbose("Identified maximum number of tracks: " + str(track_max))
        pairs_list = []
        for t in range(track_max + 1):
            ps_out_t = ps_out[ps_out["TrackID_out"] == t]
            intersect, intersect_in, intersect_out = np.intersect1d(
                ps_in[merge_columns], ps_out_t[merge_columns], return_indices=True
            )
            pairs = rfn.merge_arrays(
                (
                    intersect,
                    ps_in[ps_in_uniques][intersect_in],
                    ps_out_t[ps_out_uniques][intersect_out],
                ),
                asrecarray=True,
                flatten=True,
            )
            if len(pairs) > 0:
                pairs_list.append(pairs)
        pairs = rfn.stack_arrays(pairs_list, asrecarray=True)
        np.recarray.sort(pairs, order=["RunID", "EventID", "TrackID_in", "TrackID_out"])
    verbose("Merged input and output phase spaces.")

    if args_info.fit is not None:  # Single LUT
        verbose("Converting energy loss or TOF to WEPL using single LUT technique…")
        with open(args_info.fit, encoding="utf-8") as f:
            p = json.load(f)
        if args_info.fit_kind == "tof":
            xs = pairs["GlobalTime_out"] - pairs["GlobalTime_in"]
        elif args_info.fit_kind == "energy":
            xs = pairs["KineticEnergy_in"] - pairs["KineticEnergy_out"]
        else:
            raise NotImplementedError
        wepls = np.polyval(p, xs)
        pairs["KineticEnergy_in"] = 0.0
        pairs["KineticEnergy_out"] = wepls
    elif args_info.lut_tof is not None and args_info.lut_vel is not None:  # Double LUT
        verbose("Converting energy loss or TOF to WEPL using double LUT technique…")
        convert_tof_to_wepl(
            args_info.lut_tof,
            args_info.lut_vel,
            pairs,
            args_info.quadric,
            args_info.angle,
        )

    number_of_runs = pairs["RunID"].max() + 1
    verbose("Identified number of runs: " + str(number_of_runs))

    ComponentType = itk.ctype("float")
    PixelType = itk.Vector[ComponentType, 3]
    ImageType = itk.Image[PixelType, 2]

    run_range = range(args_info.min_run, min(number_of_runs, args_info.max_run))
    for r in run_range:
        ps_run = pairs[pairs["RunID"] == r]
        if len(ps_run) == 0:
            continue

        ps_np = np.empty(shape=(len(ps_run), 5, 3), dtype=np.float32)
        ps_np[:, 0, 0] = ps_run["u_in"]
        ps_np[:, 0, 1] = ps_run["v_in"]
        ps_np[:, 0, 2] = ps_run["w_in"]
        ps_np[:, 1, 0] = ps_run["u_out"]
        ps_np[:, 1, 1] = ps_run["v_out"]
        ps_np[:, 1, 2] = ps_run["w_out"]
        ps_np[:, 2, 0] = ps_run["du_in"]
        ps_np[:, 2, 1] = ps_run["dv_in"]
        ps_np[:, 2, 2] = ps_run["dw_in"]
        ps_np[:, 3, 0] = ps_run["du_out"]
        ps_np[:, 3, 1] = ps_run["dv_out"]
        ps_np[:, 3, 2] = ps_run["dw_out"]
        ps_np[:, 4, 0] = ps_run[measurement_column + "_in"]
        ps_np[:, 4, 1] = ps_run[measurement_column + "_out"]
        ps_np[:, 4, 2] = (
            ps_run["TrackID"] if args_info.no_nuclear else ps_run["TrackID_out"]
        )

        df_itk = itk.GetImageFromArray(ps_np, ttype=ImageType)

        output_file = args_info.output.replace(".", f"{r:04d}.")
        itk.imwrite(df_itk, output_file)
        verbose(f"Wrote file {output_file}.")


def main(argv=None):
    parser = build_parser()
    args_info = parser.parse_args(argv)
    process(args_info)


if __name__ == "__main__":
    main()
