"""Two checks of the full-range runs (0.25 to 2.55 rad/s for OC4, 0.25 to 2.50 rad/s for VolturnUS-S).

Left: the mean drift F_1 of the near field (pressure integration, MeanDriftPressure.9) against the far field (momentum,
MeanDriftMomentum.8) of the same run, |near - far| / max|far|, per frequency. The far field needs no body-surface
velocity, so the gap shows the error of the velocity method: velocity_method = source (the 0.25-1.90 rad/s OC4 run) against
velocity_method = differencing (the full-range runs, which need it because source rejects irregular-frequency removal).

Right: the median |QTF-_1| of VolturnUS-S per band of difference frequency, for WAMIT (IEA-15-240-RWT-UMaineSemi.12d,
IEA Wind Task 37, v1.1.18, Apache-2.0; free surface to RINNER = 200 m, NAL = 0), for HAMS-MREL .12d with and without its
free-surface part, and for that part alone (D_FS of QTF_DirectTerms.txt, written as .12d minus the no-FS file).

Values in N/m^2: file value x rho g with rho = 1025 kg/m^3, g = 9.80665 m/s^2 (ULEN = 1), as compare_diagonals.py.

Usage: python -I full_range_diagnostics.py OUT_PNG OUT_TXT OC4_SOURCE_DIR OC4_DIFF_DIR VOLT_DIFF_DIR IEA_12D VOLT_NOFS_12D
(the *_DIR arguments are the Wamit_format output folders)
"""
import math
import os
import statistics
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO_G = 1025.0 * 9.80665
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
BLUE, ORANGE, GREEN, PINK = "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4"


def drift(folder, name):
    out = {}
    for line in open(os.path.join(folder, name)):
        v = [float(x) for x in line.split()]
        if len(v) >= 5 and int(v[3]) == 1:
            out[round(v[0], 2)] = v[4] * RHO_G
    return out


def qtf(path, periods):
    out = {}
    for line in open(path):
        p = line.split()
        if len(p) != 9:
            continue
        v = [float(x) for x in p]
        if v[2] or v[3] or int(v[4]) != 1:
            continue
        a, b = (2 * math.pi / v[0], 2 * math.pi / v[1]) if periods else (v[0], v[1])
        a, b = round(a, 2), round(b, 2)
        c = complex(v[7], v[8]) * RHO_G
        out[(a, b)] = c
        out.setdefault((b, a), c.conjugate())
    return out


def main(out_png, out_txt, oc4_src, oc4_diff, volt_diff, iea, volt_nofs):
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.facecolor": SURFACE, "figure.facecolor": SURFACE, "legend.frameon": False,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.5, 5.2), layout="constrained")
    lines = ["Mean drift F_1: |near field - far field| / max|far field| over 0.30-1.90 rad/s and over > 1.90 rad/s, "
             "largest / median, %"]
    runs = (("OC4, 9208 panels, velocity_method = source, no lid (0.25-1.90 rad/s)", oc4_src, ORANGE),
            ("OC4, 9208 panels, velocity_method = differencing, lid (0.25-2.55 rad/s)", oc4_diff, BLUE),
            ("VolturnUS-S, 5076 panels, velocity_method = differencing, lid (0.25-2.50 rad/s)", volt_diff, GREEN))
    for label, folder, colour in runs:
        far, near = drift(folder, "MeanDriftMomentum.8"), drift(folder, "MeanDriftPressure.9")
        ws = [w for w in sorted(far) if w >= 0.30 - 1e-9]
        scale = max(abs(far[w]) for w in ws)
        err = [100 * abs(near[w] - far[w]) / scale for w in ws]
        left.plot(ws, err, color=colour, lw=2.0, label=label)
        for lo, hi in ((0.30, 1.90), (1.90 + 1e-6, 3.0)):
            e = sorted(x for w, x in zip(ws, err) if lo - 1e-9 <= w <= hi + 1e-9)
            if e:
                lines.append("  %-80s %.2f-%.2f: %6.2f / %5.2f" % (label, lo, min(hi, ws[-1]), e[-1], statistics.median(e)))
    left.set_xlabel(r"$\omega$ (rad/s)")
    left.set_ylabel(r"$|\bar F_1^{near} - \bar F_1^{far}|$ / max$|\bar F_1^{far}|$ (%)")
    left.set_title(r"Mean drift $\bar F_1$: pressure integration against momentum, same run", loc="left", color=INK)
    left.legend(fontsize=7.5, loc="upper left")

    w, h, hn = qtf(iea, True), qtf(os.path.join(volt_diff, "QTF_Direct.12d"), False), qtf(volt_nofs, False)
    edges = [0.0, 0.1, 0.2, 0.3, 0.45, 0.6, 0.75, 0.9, 1.05, 1.2, 1.5, 1.8, 2.25]
    mids, cols = [], {k: [] for k in ("wamit", "hams", "nofs", "fs")}
    lines += ["", "VolturnUS-S, median |QTF-_1| (N/m^2) per band of omega_1 - omega_2, both frequencies in 0.30-2.50 rad/s",
              "  band         WAMIT   HAMS .12d   HAMS .12d no FS   HAMS FS part   pairs"]
    for lo, hi in zip(edges[:-1], edges[1:]):
        ks = [k for k in h if k in w and k in hn and min(k) >= 0.30 - 1e-9 and max(k) <= 2.50 + 1e-9
              and lo + 1e-9 < abs(k[0] - k[1]) <= hi + 1e-9]
        if not ks:
            continue
        vals = dict(wamit=statistics.median(abs(w[k]) for k in ks), hams=statistics.median(abs(h[k]) for k in ks),
                    nofs=statistics.median(abs(hn[k]) for k in ks), fs=statistics.median(abs(h[k] - hn[k]) for k in ks))
        mids.append(0.5 * (lo + hi))
        for key in cols:
            cols[key].append(vals[key])
        lines.append("  %4.2f-%4.2f  %9.3g  %9.3g  %15.3g  %12.3g   %d" % (lo, hi, vals["wamit"], vals["hams"], vals["nofs"],
                                                                         vals["fs"], len(ks)))
    right.plot(mids, cols["wamit"], color=INK, lw=2.0, marker="o", ms=4,
               label="WAMIT v6.107S, IEA-15-240-RWT-UMaineSemi.12d (IEA Wind Task 37), free surface to 200 m")
    right.plot(mids, cols["hams"], color=BLUE, lw=2.0, marker="o", ms=4, label="HAMS-MREL .12d (direct), with free surface")
    right.plot(mids, cols["nofs"], color=PINK, lw=2.0, marker="o", ms=4, label="HAMS-MREL .12d, without the free-surface part")
    right.plot(mids, cols["fs"], color=GREEN, lw=2.0, ls=(0, (5, 2)), marker="o", ms=4,
               label="HAMS-MREL free-surface part alone (D_FS)")
    right.set_xlabel(r"$\omega_1 - \omega_2$ (rad/s), band centre")
    right.set_ylabel(r"median $|QTF^{-}_{1}|$ over the band (N/m$^2$)")
    right.set_title(r"VolturnUS-S: $QTF^{-}_{1}$ away from the diagonal", loc="left", color=INK)
    right.legend(fontsize=7.5, loc="upper left")
    fig.suptitle("Full-range runs: near-field against far-field mean drift, and the free-surface part far from the "
                 "diagonal", x=0.01, ha="left", color=INK)
    fig.savefig(out_png, dpi=150)
    open(out_txt, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(*sys.argv[1:8])
